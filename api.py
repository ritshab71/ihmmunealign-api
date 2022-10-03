import src.iHMMune_align as iHMMuneAlign

from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route('/getihmmune/', methods=['GET'])
def respond():
    sequence = request.args.get("sequence", None)

    response = {}
    if not sequence:
        response['error'] = 'No sequence found.'
    elif sequence == '':
        response['error'] = 'Empty sequence provided'
    else:
        response = iHMMuneAlign.multi_cell_align_sequence(sequence, False)

    print(response)
    return jsonify(response)

@app.route('/')
def index():
    return "<h1>Welcome to ihmmunealign-api!</h1>"

if __name__ == '__main__':
    app.run(threaded=True, port=5000)
