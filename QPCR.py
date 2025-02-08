import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def process_qpcr_data(file_path, control_samples, internal_controls, treatment_samples, output_excel, output_plot):
    # 读取 Excel 文件，将指定路径的 Excel 文件读入为 DataFrame
    df = pd.read_excel(file_path)

    # 定义一个内部函数，用于按样品分组计算每个样品的 ΔCq
    def calculate_delta_cq_per_sample(group):
        # 获取当前组对应的样品名称
        sample = group['Sample'].iloc[0]
        # 获取当前组中除内参基因外的目标基因数据
        target_data = group[~group['Target'].isin(internal_controls)]
        # 如果内参基因只有一个
        if len(internal_controls) == 1:
            # 获取当前组中内参基因的数据
            internal_data = group[group['Target'] == internal_controls[0]]
        else:
            # 如果样品在对照样品列表中，获取其索引；否则在处理样品列表中获取其索引
            index = control_samples.index(sample) if sample in control_samples else treatment_samples.index(sample)
            # 根据索引获取对应的内参基因数据
            internal_data = group[group['Target'] == internal_controls[index]]

        # 获取目标基因数据的数量
        num_target = len(target_data)
        # 获取内参基因数据的数量
        num_internal = len(internal_data)
        # 检查目标基因和内参基因的重复次数是否匹配
        if num_target % num_internal!= 0:
            # 如果不匹配，抛出异常并提示错误信息
            raise ValueError(f"样品 {sample} 中目标基因和内参基因的重复次数不匹配。")

        # 计算目标基因数据数量是内参基因数据数量的倍数
        repeats = num_target // num_internal
        # 将内参基因的 Cq 值重复相应次数，使其与目标基因数据数量匹配
        internal_cq = np.tile(internal_data['Cq'].values, repeats)
        # 获取目标基因的 Cq 值
        target_cq = target_data['Cq'].values

        # 计算目标基因与内参基因 Cq 值的差值，即 ΔCq
        delta_cq_values = target_cq - internal_cq
        # 将计算得到的 ΔCq 值添加到目标基因数据中
        target_data['ΔCq'] = delta_cq_values
        # 返回包含 ΔCq 值的目标基因数据
        return target_data

    # 按样品对数据进行分组，并调用上述函数计算每个样品的 ΔCq
    delta_cq_df = df.groupby('Sample').apply(calculate_delta_cq_per_sample).reset_index(drop=True)

    # 初始化一个列表，用于存储每个对照样品的平均 ΔCq 值 DataFrame
    control_avg_delta_cq_dfs = []
    # 遍历每个对照样品
    for control in control_samples:
        # 计算当前对照样品中每个基因的平均 ΔCq 值
        control_delta_cq_avg = delta_cq_df[delta_cq_df['Sample'] == control].groupby(['Target'])['ΔCq'].mean().reset_index()
        # 重命名平均 ΔCq 值列，使其包含对照样品名称
        control_delta_cq_avg = control_delta_cq_avg.rename(columns={'ΔCq': f'{control} Avg ΔCq'})
        # 将当前对照样品的平均 ΔCq 值 DataFrame 添加到列表中
        control_avg_delta_cq_dfs.append(control_delta_cq_avg)

    # 遍历每个对照样品的平均 ΔCq 值 DataFrame
    for control_avg_delta_cq_df in control_avg_delta_cq_dfs:
        # 将平均 ΔCq 值 DataFrame 与包含所有样品 ΔCq 值的 DataFrame 进行合并
        delta_cq_df = pd.merge(delta_cq_df, control_avg_delta_cq_df, on='Target', how='left')

    # 定义一个函数，用于计算每个样品的 ΔΔCq 值
    def calculate_delta_delta_cq(row):
        # 遍历每个对照样品
        for control in control_samples:
            # 如果当前样品名称中包含对照样品的时间信息
            if control.split('-')[-1] in row['Sample']:
                # 计算当前样品的 ΔΔCq 值
                return row['ΔCq'] - row[f'{control} Avg ΔCq']
        # 如果未匹配到对照样品，返回 NaN
        return np.nan

    # 对包含所有样品 ΔCq 值的 DataFrame 应用上述函数，计算每个样品的 ΔΔCq 值
    delta_cq_df['ΔΔCq'] = delta_cq_df.apply(calculate_delta_delta_cq, axis=1)

    # 计算每个样品的相对表达量，通过公式 2 ** -ΔΔCq
    delta_cq_df['Relative Expression'] = 2 ** -(delta_cq_df['ΔΔCq'])

    # 初始化一个列表，用于存储每个样品的标准差信息
    stdev_list = []
    # 合并对照样品和处理样品列表
    all_samples = control_samples + treatment_samples
    # 遍历每个目标基因
    for target in delta_cq_df['Target'].unique():
        # 如果目标基因不是内参基因
        if target not in internal_controls:
            # 遍历每个样品
            for sample in all_samples:
                # 获取当前目标基因和样品对应的相对表达量数据
                fold_changes = delta_cq_df[(delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == sample)][
                    'Relative Expression'].values
                # 计算相对表达量数据的标准差
                stdev = np.std(fold_changes)
                # 将当前目标基因、样品和标准差信息添加到列表中
                stdev_list.append({
                    'Target': target,
                    'Sample': sample,
                    'STDEV': stdev
                })

    # 将标准差信息列表转换为 DataFrame
    stdev_df = pd.DataFrame(stdev_list)

    # 按目标基因和样品对相对表达量数据进行分组，并计算平均值
    result_df = delta_cq_df.groupby(['Target', 'Sample']).agg({'Relative Expression': 'mean'}).reset_index()
    # 将平均值 DataFrame 与标准差 DataFrame 进行合并
    result_df = pd.merge(result_df, stdev_df, on=['Target', 'Sample'])

    # 从样品名称中提取时间信息（如 2h 或 4h）
    result_df['Time'] = result_df['Sample'].str.extract(r'(\d+h)')
    # 按目标基因、时间和样品对结果 DataFrame 进行排序
    result_df = result_df.sort_values(by=['Target', 'Time', 'Sample'])
    # 删除时间信息列
    result_df = result_df.drop(columns=['Time'])

    # 初始化一个列表，用于存储 T - test 的结果
    t_test_results = []
    # 遍历每个目标基因
    for target in delta_cq_df['Target'].unique():
        # 如果目标基因不是内参基因
        if target not in internal_controls:
            # 遍历对照样品和处理样品，一一对应
            for control, treatment in zip(control_samples, treatment_samples):
                # 获取当前目标基因和对照样品对应的相对表达量数据
                control_relative_expression = delta_cq_df[
                    (delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == control)]['Relative Expression'].values
                # 获取当前目标基因和处理样品对应的相对表达量数据
                treatment_relative_expression = delta_cq_df[
                    (delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == treatment)]['Relative Expression'].values
                # 进行独立样本 T - test，计算 T 统计量和 p 值
                t_stat, p_value = stats.ttest_ind(control_relative_expression, treatment_relative_expression)
                # 将当前目标基因、时间点、T 统计量和 p 值添加到结果列表中
                t_test_results.append({
                    'Target': target,
                    'Time Point': control.split('-')[-1],
                    'T - test Statistic': t_stat,
                    'P - value': p_value
                })

    # 将 T - test 结果列表转换为 DataFrame
    t_test_df = pd.DataFrame(t_test_results)

    # 使用 ExcelWriter 将结果保存到 Excel 文件中
    with pd.ExcelWriter(output_excel) as writer:
        # 将相对表达量结果写入 Excel 文件的 'Relative Expression' 工作表
        result_df.to_excel(writer, sheet_name='Relative Expression', index=False)
        # 将 T - test 结果写入 Excel 文件的 'T - test Results' 工作表
        t_test_df.to_excel(writer, sheet_name='T - test Results', index=False)

    # 打印提示信息，表明处理完成并保存到 Excel 文件
    print(f"处理完成，结果已保存到 {output_excel} 文件中。")

    # 创建一个新的图形，设置图形大小
    plt.figure(figsize=(12, 8))
    # 使用 seaborn 的 barplot 绘制相对表达量的柱状图，不显示置信区间
    bars = sns.barplot(x='Target', y='Relative Expression', hue='Sample', data=result_df, errorbar=None)

    # 获取结果 DataFrame 中目标基因的数量
    num_targets = len(result_df['Target'].unique())
    # 获取结果 DataFrame 中样品的数量
    num_samples = len(result_df['Sample'].unique())
    # 计算柱状图中每个柱子的宽度
    width = 0.8 / num_samples
    # 遍历每个目标基因
    for i, target in enumerate(result_df['Target'].unique()):
        # 遍历每个样品
        for j, sample in enumerate(result_df['Sample'].unique()):
            # 获取当前目标基因和样品对应的结果数据子集
            subset = result_df[(result_df['Target'] == target) & (result_df['Sample'] == sample)]
            # 如果子集不为空
            if not subset.empty:
                # 获取当前样品的相对表达量平均值
                rel_expr = subset['Relative Expression'].values[0]
                # 获取当前样品的相对表达量标准差
                stdev = subset['STDEV'].values[0]
                # 调整误差线的 x 位置，使其与柱子中心对齐
                x_pos = i - 0.4 + (j + 0.5) * width
                # 绘制误差线
                plt.errorbar(x_pos, rel_expr, yerr=stdev, fmt='none', color='black', capsize=5)

    # 遍历 T - test 结果 DataFrame 的每一行
    for _, row in t_test_df.iterrows():
        # 获取当前行的目标基因
        target = row['Target']
        # 获取当前行的时间点
        time_point = row['Time Point']
        # 获取当前行的 p 值
        p_value = row['P - value']
        # 遍历对照样品和处理样品，一一对应
        for control, treatment in zip(control_samples, treatment_samples):
            # 如果时间点在对照样品名称中
            if time_point in control:
                # 设置当前处理样品
                sample = treatment
                # 设置当前对照样品
                control_sample = control
                # 跳出循环
                break
        # 获取目标基因在结果 DataFrame 中目标基因列表的索引
        target_index = list(result_df['Target'].unique()).index(target)
        # 获取处理样品在结果 DataFrame 中样品列表的索引
        sample_index = list(result_df['Sample'].unique()).index(sample)
        # 获取当前目标基因和处理样品对应的结果数据子集
        subset = result_df[(result_df['Target'] == target) & (result_df['Sample'] == sample)]
        # 如果子集不为空
        if not subset.empty:
            # 获取当前样品的相对表达量平均值
            rel_expr = subset['Relative Expression'].values[0]
            # 获取当前样品的相对表达量标准差
            stdev = subset['STDEV'].values[0]
            # 计算显著性标记的 x 位置
            x_pos = target_index - 0.4 + (sample_index + 0.5) * width
            # 计算显著性标记的 y 位置，设置在误差线顶端
            y_pos = rel_expr + stdev
            # 如果 p 值小于 0.01
            if p_value < 0.01:
                # 在误差线顶端绘制 '**' 标记，表示显著性
                plt.text(x_pos, y_pos + 0.1, '**', ha='center', va='bottom')
            # 如果 p 值小于 0.05 且大于等于 0.01
            elif p_value < 0.05:
                # 在误差线顶端绘制 '*' 标记，表示显著性
                plt.text(x_pos, y_pos + 0.1, '*', ha='center', va='bottom')

    # 遍历每个目标基因
    for i, target in enumerate(result_df['Target'].unique()):
        # 遍历每个样品
        for j, sample in enumerate(result_df['Sample'].unique()):
            # 获取当前目标基因和样品对应的原始相对表达量数据
            individual_values = delta_cq_df[(delta_cq_df['Target'] == target) & (delta_cq_df['Sample'] == sample)][
                'Relative Expression'].values
            # 计算绘制点的 x 位置
            x_pos = i - 0.4 + (j + 0.5) * width
            # 遍历每个原始相对表达量数据点
            for value in individual_values:
                # 绘制空心黑色圆点
                plt.plot(x_pos, value, 'ko', fillstyle='none', markersize=3)

    # 设置 x 轴标签
    plt.xlabel('Genes')
    # 设置 y 轴标签
    plt.ylabel('Relative Expression')
    # 设置图形标题
    plt.title('Relative Expression of Genes in Different Samples')
    # 设置图例标题
    plt.legend(title='Samples')
    # 自动调整子图参数，使其填充整个图形区域
    plt.tight_layout()
    # 保存图形为指定路径的图片文件
    plt.savefig(output_plot)
    # 打印提示信息，表明绘图完成并保存到图片文件
    print(f"绘图完成，结果已保存到 {output_plot} 文件中。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='处理 qPCR 数据')
    parser.add_argument('input_file', help='输入的 Excel 文件路径')
    parser.add_argument('-c', '--control_samples', nargs='+', required=True, help='对照样品名称')
    parser.add_argument('-i', '--internal_controls', nargs='+', required=True, help='内参基因名称')
    parser.add_argument('-t', '--treatment_samples', nargs='+', required=True, help='处理样品名称')
    # 设置默认输出文件名
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
        # 如果只提供了一个值，设置默认的图片文件名
        output_excel = args.output_files
        output_plot ='result.png'

    process_qpcr_data(input_file, control_samples, internal_controls, treatment_samples, output_excel, output_plot)
