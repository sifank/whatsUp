import sys

# Add your application directory to the system path
sys.path.insert(0, "/var/www/whatsUp")

# Import the 'app' object from your app.py file
from whatsUp import app as application

