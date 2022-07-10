import src.iHMMune_align as iHMMuneAlign

from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route('/getihmmune/', methods=['GET'])
def respond():
    input_file = 'src/files/input.txt'
    output = iHMMuneAlign.multi_cell_align(input_file, False)
    return jsonify(output)

@app.route('/')
def index():
    return "<h1>Welcome to ihmmunealign-api!</h1>"

if __name__ == '__main__':
    app.run(threaded=True, port=5000)
