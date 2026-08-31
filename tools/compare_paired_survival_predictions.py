import argparse
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from dataset.label_converter import MetaSurvData
from eval.SurvivalEVAL.Evaluations.Concordance import concordance
from eval.SurvivalEVAL.Evaluations.util import predict_mean_survival_time
from utils.io import read_file_data_splitting


def resolve_project_path(path):
    path = Path(path)
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_predictions(path):
    frame = pd.read_csv(path, dtype={'patient_id': str})
    survival_columns = [column for column in frame.columns if column.startswith('surf_')]
    required = {'patient_id', 'e', *survival_columns}
    missing = required - set(frame.columns)
    if missing or not survival_columns:
        raise ValueError(f'{path} is missing patient IDs, event labels, or survival curves.')
    if frame['patient_id'].duplicated().any():
        raise ValueError(f'{path} contains duplicate patient IDs.')
    survival = frame[survival_columns].to_numpy(dtype=float)
    if not np.isfinite(survival).all():
        raise ValueError(f'{path} contains non-finite survival probabilities.')
    frame = frame[['patient_id', 'e']].copy()
    frame['row_index'] = np.arange(len(frame))
    return frame, survival, survival_columns


def load_config(prediction_path, explicit_path=None):
    config_path = Path(explicit_path) if explicit_path else Path(prediction_path).parent / 'config.yaml'
    if not config_path.exists():
        raise FileNotFoundError(f'Cannot find result config: {config_path}')
    with config_path.open('r', encoding='utf-8') as stream:
        return yaml.safe_load(stream), config_path


def build_meta(config):
    split_path = str(config['data_split_path'])
    if '{' in split_path:
        raise ValueError(
            'The saved config still contains data-split placeholders. '
            'Pass a resolved result config with --config.'
        )
    data_split = read_file_data_splitting(str(resolve_project_path(split_path)))
    meta = MetaSurvData(
        str(resolve_project_path(config['path_table'])),
        data_split=data_split,
        verbose=False,
    )
    if config.get('time_format', 'interval') == 'interval':
        meta.generate_discrete_label(
            num_bins=config.get('time_bins'),
            new_column_t='y_t',
            new_column_e='y_e',
            use_quantiles=False,
            summary=False,
            open_ended_last_bin=bool(config.get('open_ended_last_time_bin', False)),
        )
    else:
        raise ValueError('This comparison currently supports interval survival outputs only.')
    return meta


def predicted_event_times(survival, time_coordinates):
    return np.asarray([
        predict_mean_survival_time(curve, time_coordinates, interpolation='Linear')
        for curve in survival
    ], dtype=float)


def c_index(predicted_times, event_times, event_indicators, train_times, train_indicators):
    return float(concordance(
        predicted_times,
        event_times,
        event_indicators,
        train_times,
        train_indicators,
        pair_method='Comparable',
        ties='All',
    )[0])


def prediction_split(path):
    match = re.search(r'_pred_(train|validation|test)\.csv$', Path(path).name)
    return match.group(1) if match else None


def read_logged_metric(prediction_path, split):
    if split is None:
        return None
    metrics_path = Path(prediction_path).parent / 'train_metrics-last.txt'
    if not metrics_path.exists():
        return None
    pattern = re.compile(rf'^{re.escape(split)}/pred_c_index\s+-->\s+([-+0-9.eE]+)')
    for line in metrics_path.read_text(encoding='utf-8').splitlines():
        match = pattern.search(line.strip())
        if match:
            return float(match.group(1))
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Paired bootstrap using the exact continuous-time SurvivalEVAL C-index contract.'
    )
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--candidate', required=True)
    parser.add_argument('--config', default=None)
    parser.add_argument('--bootstrap', type=int, default=5000)
    parser.add_argument('--seed', type=int, default=2404)
    parser.add_argument('--metric_tolerance', type=float, default=1e-6)
    parser.add_argument('--output', default=None)
    args = parser.parse_args()

    baseline, baseline_survival, baseline_columns = load_predictions(args.baseline)
    candidate, candidate_survival, candidate_columns = load_predictions(args.candidate)
    if baseline_columns != candidate_columns:
        raise ValueError('Baseline and candidate survival-curve columns differ.')

    paired = baseline.merge(candidate, on='patient_id', how='inner', suffixes=('_base', '_cand'))
    if len(paired) != len(baseline) or len(paired) != len(candidate):
        raise ValueError('Baseline and candidate patient sets are not identical.')
    if not np.allclose(paired['e_base'], paired['e_cand']):
        raise ValueError('Baseline and candidate event labels differ.')

    config, config_path = load_config(args.baseline, args.config)
    candidate_config, candidate_config_path = load_config(args.candidate)
    contract_keys = (
        'path_table', 'data_split_path', 'time_bins', 'time_format',
        'open_ended_last_time_bin',
    )
    for key in contract_keys:
        if config.get(key) != candidate_config.get(key):
            raise ValueError(
                f'Evaluation setting {key!r} differs between {config_path} and '
                f'{candidate_config_path}.'
            )

    meta = build_meta(config)
    actual = meta.get_patient_data(
        pids=paired['patient_id'].tolist(), ret_columns=['t', 'e']
    )
    train = meta.get_patient_data(split='train', ret_columns=['t', 'e'])
    if not np.allclose(actual.e.to_numpy(dtype=float), paired['e_base'].to_numpy(dtype=float)):
        raise ValueError('Saved event labels do not match the continuous-time metadata.')

    baseline_times_all = predicted_event_times(baseline_survival, meta.time_coordinates)
    candidate_times_all = predicted_event_times(candidate_survival, meta.time_coordinates)
    baseline_times = baseline_times_all[paired['row_index_base'].to_numpy(dtype=int)]
    candidate_times = candidate_times_all[paired['row_index_cand'].to_numpy(dtype=int)]
    event_times = actual.t.to_numpy(dtype=float)
    event_indicators = actual.e.to_numpy(dtype=float)
    train_times = train.t.to_numpy(dtype=float)
    train_indicators = train.e.to_numpy(dtype=float)

    baseline_c = c_index(
        baseline_times, event_times, event_indicators, train_times, train_indicators
    )
    candidate_c = c_index(
        candidate_times, event_times, event_indicators, train_times, train_indicators
    )

    split = prediction_split(args.baseline)
    if split != prediction_split(args.candidate):
        raise ValueError('Baseline and candidate prediction splits differ.')
    logged_baseline = read_logged_metric(args.baseline, split)
    logged_candidate = read_logged_metric(args.candidate, split)
    for name, computed, logged in (
        ('baseline', baseline_c, logged_baseline),
        ('candidate', candidate_c, logged_candidate),
    ):
        if logged is not None and abs(computed - logged) > args.metric_tolerance:
            raise RuntimeError(
                f'{name} C-index contract mismatch: recomputed={computed:.12f}, '
                f'logged={logged:.12f}. Refusing to bootstrap incomparable metrics.'
            )

    observed_delta = candidate_c - baseline_c
    rng = np.random.RandomState(args.seed)
    deltas = []
    size = len(paired)
    for _ in range(args.bootstrap):
        sampled = rng.randint(0, size, size=size)
        sampled_events = event_indicators[sampled]
        if not np.any(sampled_events):
            continue
        try:
            deltas.append(
                c_index(
                    candidate_times[sampled], event_times[sampled], sampled_events,
                    train_times, train_indicators,
                )
                - c_index(
                    baseline_times[sampled], event_times[sampled], sampled_events,
                    train_times, train_indicators,
                )
            )
        except (ValueError, ZeroDivisionError):
            continue
    if not deltas:
        raise RuntimeError('No valid bootstrap samples contained comparable pairs.')

    deltas = np.asarray(deltas, dtype=float)
    summary = pd.DataFrame([{
        'metric_contract': 'SurvivalEVAL_mean_continuous_comparable_all_ties',
        'n': size,
        'events': int(event_indicators.sum()),
        'baseline_c_index': baseline_c,
        'candidate_c_index': candidate_c,
        'observed_delta': observed_delta,
        'bootstrap_samples': len(deltas),
        'delta_ci025': float(np.quantile(deltas, 0.025)),
        'delta_ci975': float(np.quantile(deltas, 0.975)),
        'probability_delta_gt_zero': float(np.mean(deltas > 0.0)),
        'logged_baseline_c_index': logged_baseline,
        'logged_candidate_c_index': logged_candidate,
    }])
    print(summary.to_string(index=False, float_format=lambda value: f'{value:.6f}'))

    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(output, index=False)
        print(f'output={output}')


if __name__ == '__main__':
    main()
