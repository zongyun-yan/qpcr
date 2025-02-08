from flask import Flask, request, render_template, send_file
import pandas as pd
import os
import QPCR
import tempfile

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/process', methods=['POST'])
def process():
    # 获取上传文件
    file = request.files['file']
    # 创建临时文件存储上传数据
    with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as temp:
        file.save(temp.name)
        temp.close()
    # 获取表单参数
    control_samples = request.form.get('control_samples').split(',')
    internal_controls = request.form.get('internal_controls').split(',')
    treatment_samples = request.form.get('treatment_samples').split(',')
    # 生成输出文件路径
    output_excel = 'output.xlsx'
    output_plot = 'output.png'
    try:
        # 调用qPCR数据处理函数
        QPCR.process_qpcr_data(temp.name, control_samples, internal_controls, treatment_samples, output_excel, output_plot)
        # 返回处理结果文件
        return send_file(output_excel, as_attachment=True)
    except Exception as e:
        return f"处理出错: {str(e)}"
    finally:
        # 处理完成后删除临时文件
        os.remove(temp.name)


if __name__ == '__main__':
    app.run(debug=True)