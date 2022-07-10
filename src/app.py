import iHMMune_align as iHMMuneAlign
import mutability_score as MutabilityScore

from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route('/getihmmune/', methods=['GET'])
def respond():
    input_file = 'files/input.txt'
    output = iHMMuneAlign.multi_cell_align(input_file, False)
    return jsonify(output)

if __name__ == '__main__':
    app.run(threaded=True, port=5000)
