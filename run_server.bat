@echo off
echo Running Flask server...
git pull origin main
python app/router.py
pause
