from flask import Flask, render_template, request, redirect
import csv
import os

app = Flask(__name__)

LOG_PATH = os.path.join(os.path.dirname(__file__), '..', 'database', 'logs.csv')

@app.route('/')
def index():
    return render_template('sales.html')

@app.route('/submit', methods=['POST'])
def submit():
    name = request.form['name']
    order_number = request.form['order_number']
    item = request.form['item']
    price = request.form['price']
    quantity = request.form['quantity']
    total_amount = request.form['total_amount']

    with open(LOG_PATH, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([name, order_number, item, price, quantity, total_amount])

    return redirect('/')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
