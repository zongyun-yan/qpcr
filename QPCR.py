import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os


def process_qpcr_data(file_path, control_samples, internal_controls, treatment_samples, output_excel, output_plot):
    df = pd.read_excel(file_path)
    required_columns = ['Sample', 'Target', 'Cq']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"输入的Excel文件缺少必要的列: {', '.join(missing_columns)}")

    def calculate_delta_cq_per_sample(group):
        sample = group['Sample'].iloc[0]
        target_data = group[~group['Target'].isin(internal_controls)]
        if len(internal_controls) == 1:
            internal_data = group[group['Target'] == internal_controls[0]]
        else:
            index = control_samples.index(sample) if sample in control_samples else treatment_samples.index(sample)
            internal_data = group[group['Target'] == internal_controls[index]]
        num_target = len(target_data)
        num_internal = len(internal_data)
        if num_target % num_internal!= 0:
            raise ValueError(f"样品 {sample} 中目标基因和内参基因的重复次数不匹配。")
        repeats = num_target // num_internal
        internal_cq = np.tile(internal_data['Cq'].values, repeats)
        target_cq = target_data['Cq'].values
        delta_cq_values = target_cq - internal_cq
        target_data['ΔCq'] = delta_cq_values
        return target_data

    delta_cq_df = df.groupby('Sample').apply(calculate_delta_cq_per_sample).reset_index(drop=True)
    control_avg_delta_cq_dfs = []
    for control in control_samples:
        control_delta_cq_avg = delta_cq_df[delta_cq_df['Sample'] == control].groupby(['Target'])['ΔCq'].mean().reset_index()
        control_delta_cq_avg = control_delta_cq_avg.rename(columns={'ΔCq': f'{control} Avg ΔCq'})
        control_avg_delta_cq_dfs.append(control_delta_cq_avg)
    for control_avg_delta_cq_df in control_avg_delta_cq_dfs:
        delta_cq_df = pd.merge(delta_cq_df, control_avg_delta_cq_df, on='Target', how='left')

    def calculate_delta_delta_cq(row):
        for control in control_samples:
            if control.split('-')[-1] in row['Sample']:
                return row['ΔCq'] - row[f'{control} Avg ΔCq']
        return np.nan

    delta_cq_df['ΔΔCq'] = delta_cq_df.apply(calculate_delta_delta_cq, axis=1)
    delta_cq_df['Relative Expression'] = 2 ** -(delta_cq_df['ΔΔCq'])
    stdev_list = []
    all_samples = control_samples + treatment_samples
    for target in delta_cq_df['Target'].unique():
        if target not in internal_controls:
            for sample in all_samples:
                fold_changes = delta_cq_df[(delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == sample)][
                    'Relative Expression'].values
                stdev = np.std(fold_changes)
                stdev_list.append({
                    'Target': target,
                    'Sample': sample,
                    'STDEV': stdev
                })
    stdev_df = pd.DataFrame(stdev_list)
    result_df = delta_cq_df.groupby(['Target', 'Sample']).agg({'Relative Expression':'mean'}).reset_index()
    result_df = pd.merge(result_df, stdev_df, on=['Target', 'Sample'])
    result_df['Time'] = result_df['Sample'].str.extract(r'(\d+h)')
    result_df = result_df.sort_values(by=['Target', 'Time', 'Sample'])
    result_df = result_df.drop(columns=['Time'])

    t_test_results = []
    for target in delta_cq_df['Target'].unique():
        if target not in internal_controls:
            for control, treatment in zip(control_samples, treatment_samples):
                control_relative_expression = delta_cq_df[
                    (delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == control)]['Relative Expression'].values
                treatment_relative_expression = delta_cq_df[
                    (delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == treatment)]['Relative Expression'].values
                t_stat, p_value = stats.ttest_ind(control_relative_expression, treatment_relative_expression)
                t_test_results.append({
                    'Target': target,
                    'Time Point': control.split('-')[-1],
                    'T - test Statistic': t_stat,
                    'P - value': p_value
                })
    t_test_df = pd.DataFrame(t_test_results)

    with pd.ExcelWriter(output_excel) as writer:
        result_df.to_excel(writer, sheet_name='Relative Expression', index=False)
        t_test_df.to_excel(writer, sheet_name='T - test Results', index=False)
    print(f"处理完成，结果已保存到 {output_excel} 文件中。")

    plt.figure(figsize=(12, 8))
    bars = sns.barplot(x='Target', y='Relative Expression', hue='Sample', data=result_df, errorbar=None)
    num_targets = len(result_df['Target'].unique())
    num_samples = len(result_df['Sample'].unique())
    width = 0.8 / num_samples
    for i, target in enumerate(result_df['Target'].unique()):
        for j, sample in enumerate(result_df['Sample'].unique()):
            subset = result_df[(result_df['Target'] == target) & (result_df['Sample'] == sample)]
            if not subset.empty:
                rel_expr = subset['Relative Expression'].values[0]
                stdev = subset['STDEV'].values[0]
                x_pos = i - 0.4 + (j + 0.5) * width
                plt.errorbar(x_pos, rel_expr, yerr=stdev, fmt='none', color='black', capsize=5)
    for _, row in t_test_df.iterrows():
        target = row['Target']
        time_point = row['Time Point']
        p_value = row['P - value']
        for control, treatment in zip(control_samples, treatment_samples):
            if time_point in control:
                sample = treatment
                control_sample = control
                break
        target_index = list(result_df['Target'].unique()).index(target)
        sample_index = list(result_df['Sample'].unique()).index(sample)
        subset = result_df[(result_df['Target'] == target) & (result_df['Sample'] == sample)]
        if not subset.empty:
            rel_expr = subset['Relative Expression'].values[0]
            stdev = subset['STDEV'].values[0]
            x_pos = target_index - 0.4 + (sample_index + 0.5) * width
            y_pos = rel_expr + stdev
            if p_value < 0.01:
                plt.text(x_pos, y_pos + 0.1, '**', ha='center', va='bottom')
            elif p_value < 0.05:
                plt.text(x_pos, y_pos + 0.1, '*', ha='center', va='bottom')
    for i, target in enumerate(result_df['Target'].unique()):
        for j, sample in enumerate(result_df['Sample'].unique()):
            individual_values = delta_cq_df[(delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == sample)][
                'Relative Expression'].values
            x_pos = i - 0.4 + (j + 0.5) * width
            for value in individual_values:
                plt.plot(x_pos, value, 'ko', fillstyle='none', markersize=3)
    plt.xlabel('Genes')
    plt.ylabel('Relative Expression')
    plt.title('Relative Expression of Genes in Different Samples')
    plt.legend(title='Samples')
    plt.tight_layout()
    plt.savefig(output_plot)
    print(f"绘图完成，结果已保存到 {output_plot} 文件中。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='处理 qPCR 数据')
    parser.add_argument('input_file', help='输入的 Excel 文件路径')
    parser.add_argument('-c', '--control_samples', nargs='+', required=True, help='对照样品名称')
    parser.add_argument('-i', '--internal_controls', nargs='+', required=True, help='内参基因名称')
    parser.add_argument('-t', '--treatment_samples', nargs='+', required=True, help='处理样品名称')
    parser.add_argument('-o', '--output_files', nargs='?', default=['result.xlsx','result.png'],
                        help='输出的 Excel 和图片文件路径')
    args = parser.parse_args()
    input_file = args.input_file
    control_samples = args.control_samples
    internal_controls = args.internal_controls
    treatment_samples = args.treatment_samples
    if isinstance(args.output_files, list):
        output_excel, output_plot = args.output_files
    else:
        output_excel = args.output_files
        output_plot ='result.png'
    process_qpcr_data(input_file, control_samples, internal_controls, treatment_samples, output_excel, output_plot)
