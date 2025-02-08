import os
from flask import Flask, request, render_template, send_file
import pandas as pd
import tempfile
import QPCR

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/process', methods=['POST'])
def process():
    file = request.files['file']
    with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as temp:
        file.save(temp.name)
        temp.close()
    control_samples = [sample.strip() for sample in request.form.get('control_samples').split(',')]
    internal_controls = [control.strip() for control in request.form.get('internal_controls').split(',')]
    treatment_samples = [sample.strip() for sample in request.form.get('treatment_samples').split(',')]
    current_dir = os.path.dirname(os.path.abspath(__file__))
    output_excel = os.path.join(current_dir, 'output.xlsx')
    output_plot = os.path.join(current_dir, 'output.png')
    try:
        if not os.access(current_dir, os.W_OK):
            raise PermissionError(errno.EACCES, os.strerror(errno.EACCES), current_dir)
        QPCR.process_qpcr_data(temp.name, control_samples, internal_controls, treatment_samples, output_excel, output_plot)
        return send_file(output_excel, as_attachment=True)
    except PermissionError as pe:
        return f"权限不足: {str(pe)}"
    except Exception as e:
        return f"处理出错: {str(e)}"
    finally:
        os.remove(temp.name)


@app.route('/download_plot', methods=['GET'])
def download_plot():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    output_plot = os.path.join(current_dir, 'output.png')
    try:
        return send_file(output_plot, as_attachment=True)
    except FileNotFoundError:
        return f"图片文件未找到，可能处理过程中出现异常，未成功生成图片。", 404
    except Exception as e:
        return f"下载图片出错: {str(e)}", 500


if __name__ == '__main__':
    app.run(debug=True)
