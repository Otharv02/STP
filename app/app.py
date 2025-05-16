
from flask import Flask, render_template, request, redirect, url_for, jsonify
import csv
import os
from datetime import datetime

app = Flask(__name__)

# Ensure the database directory exists
os.makedirs('database', exist_ok=True)
LOG_FILE = 'logs.csv'  # Changed to use the logs.csv file directly

def write_log(customer_info, items):
    file_exists = os.path.isfile(LOG_FILE)

    with open(LOG_FILE, mode='a', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)

        if not file_exists:
            # Create header row with customer info fields
            header_row = [
                'Name', 'Phone', 'Order By', 'Invoice Number', 'Invoice Date',
                'Due Date', 'Reference', 'Round Off', 'Total'
            ]
            
            # Add item headers for each of the 10 possible items
            for i in range(1, 11):
                header_row.extend([
                    f'Item{i}_Code', f'Item{i}_Category', f'Item{i}_Description',
                    f'Item{i}_Quantity', f'Item{i}_Unit', f'Item{i}_Price',
                    f'Item{i}_Amount'
                ])
                
            writer.writerow(header_row)

        row = [
            customer_info.get('name', ''),
            customer_info.get('phone', ''),
            customer_info.get('order_by', ''),
            customer_info.get('invoice_number', ''),
            customer_info.get('invoice_date', ''),
            customer_info.get('due_date', ''),
            customer_info.get('reference', ''),
            customer_info.get('round_off', ''),
            customer_info.get('total', '')
        ]

        for i in range(10):
            if i < len(items):
                item = items[i]
                row.extend([
                    item.get('item_code', ''),
                    item.get('category', ''),
                    item.get('description', ''),
                    item.get('quantity', ''),
                    item.get('unit', ''),
                    item.get('price', ''),
                    item.get('amount', '')
                ])
            else:
                row.extend([''] * 7)

        writer.writerow(row)

def read_logs(name_filter=None):
    if not os.path.isfile(LOG_FILE):
        return []

    logs = []
    try:
        with open(LOG_FILE, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.reader(file)
            # Skip the header row
            try:
                next(reader)
            except StopIteration:
                return []  # Return empty list if file is empty
            
            for row in reader:
                if not row or len(row) < 9:  # Make sure there are minimum required fields
                    continue
                
                # Extract basic customer info
                name = row[0]
                phone = row[1]
                order_by = row[2]
                invoice_number = row[3]
                invoice_date = row[4]
                due_date = row[5]
                reference = row[6]
                round_off = row[7]
                total = row[8]
                
                # Apply name filter if provided
                if name_filter and name_filter.lower() not in name.lower():
                    continue
                
                # Process items (make sure we're accessing valid indices)
                items_list = []  # Changed from 'items' to 'items_list' to avoid naming conflict
                for i in range(10):  # Process up to 10 items
                    base_idx = 9 + (i * 7)  # Start of item fields for item i
                    
                    if base_idx < len(row) and row[base_idx]:  # Check if item code exists
                        item = {
                            'code': row[base_idx],
                            'category': row[base_idx + 1] if base_idx + 1 < len(row) else '',
                            'description': row[base_idx + 2] if base_idx + 2 < len(row) else '',
                            'quantity': row[base_idx + 3] if base_idx + 3 < len(row) else '',
                            'unit': row[base_idx + 4] if base_idx + 4 < len(row) else '',
                            'price': row[base_idx + 5] if base_idx + 5 < len(row) else '',
                            'amount': row[base_idx + 6] if base_idx + 6 < len(row) else ''
                        }
                        items_list.append(item)
                
                # Create the log entry
                log_entry = {
                    'Name': name,
                    'Phone': phone,
                    'Order By': order_by,
                    'Invoice Number': invoice_number,
                    'Invoice Date': invoice_date,
                    'Due Date': due_date,
                    'Reference': reference,
                    'Round Off': round_off,
                    'Total': total,
                    'items_list': items_list  # Changed from 'items' to 'items_list'
                }
                
                logs.append(log_entry)
    except Exception as e:
        print(f"Error reading logs: {e}")
        import traceback
        traceback.print_exc()
    
    return logs

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/sales')
def sales_form():
    timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
    invoice_number = f"INV-{timestamp}"
    return render_template('sales.html', invoice_number=invoice_number)

@app.route('/submit', methods=['POST'])
def submit():
    customer_info = {
        'name': request.form.get('name', ''),
        'phone': request.form.get('phone', ''),
        'order_by': request.form.get('order_by', ''),
        'invoice_number': request.form.get('invoice_number', ''),
        'invoice_date': request.form.get('invoice_date', ''),
        'due_date': request.form.get('due_date', ''),
        'reference': request.form.get('reference', ''),
        'round_off': request.form.get('round_off', '0'),
        'total': request.form.get('total', '0')
    }

    total_rows = int(request.form.get('total_rows', 1))
    items = []
    for i in range(1, total_rows + 1):
        item_code = request.form.get(f'item_code_{i}', '')
        if item_code:
            items.append({
                'item_code': item_code,
                'category': request.form.get(f'category_{i}', ''),
                'description': request.form.get(f'description_{i}', ''),
                'quantity': request.form.get(f'quantity_{i}', ''),
                'unit': request.form.get(f'unit_{i}', ''),
                'price': request.form.get(f'price_{i}', ''),
                'amount': request.form.get(f'amount_{i}', '')
            })

    write_log(customer_info, items)
    return redirect(url_for('sales_form'))

@app.route('/dashboard')
def dashboard():
    try:
        name_filter = request.args.get('name', '')
        logs = read_logs(name_filter)
        return render_template('dashboard.html', logs=logs, name_filter=name_filter)
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"Error in dashboard route: {e}")
        print(error_details)
        return f"An error occurred: {str(e)}<br><pre>{error_details}</pre>", 500

@app.route('/api/search', methods=['GET'])
def search_api():
    try:
        name_filter = request.args.get('name', '')
        logs = read_logs(name_filter)
        return jsonify(logs)
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"Error in search API: {e}")
        print(error_details)
        return jsonify({"error": str(e), "details": error_details}), 500

if __name__ == '__main__':
    app.run(debug=True)