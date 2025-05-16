@echo off
echo Pulling latest code from Git...
git pull

echo Starting Flask server...
REM Set the Flask app environment variable (optional)
set FLASK_APP=app.py

REM Run the Flask app
python app\app.py

pause
