#!/usr/bin/python3

from flask import Flask, render_template, request
from datetime import datetime

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('enterDate.html')

@app.route('/submit-dropdown', methods=['POST'])
def submit_dropdown():
    # 1. Collect calendar string and pull-down strings
    date_str = request.form.get('user_date')      # e.g., "2026-07-10"
    hour_str = request.form.get('user_hour')      # e.g., "03"
    minute_str = request.form.get('user_minute')  # e.g., "30"
    period_str = request.form.get('user_period')  # e.g., "PM"

    # 2. Combine individual time fragments into a single string
    time_str = f"{hour_str}:{minute_str} {period_str}" # e.g., "03:30 PM"

    # 3. Parse individual components into date and time objects
    parsed_date = datetime.strptime(date_str, '%Y-%m-%d').date()
    parsed_time = datetime.strptime(time_str, '%I:%M %p').time() # %I is for 12-hour clock

    # 4. Merge into a final datetime object
    final_datetime = datetime.combine(parsed_date, parsed_time)

    formatted_result = final_datetime.strftime('%A, %B %d, %Y at %I:%M %p')
    return f"Drop-down choice saved for: {formatted_result}"


if __name__ == '__main__':
    app.run(host='0.0.0.0',port=5011, debug=True)

