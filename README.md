# whatsUp
Comprehensive tool for scheduling nighttime observations

cd /var/www/html
Run: git clone https://github.com/sifank/whatsUp.git
Move whatsUp.conf to /etc/apache2/sites-available
chown -R www-data:www-data /var/www/whatsUp
add "Listen 5004" to /etc/apache2/ports.conf
sudo ufw allow 5004/tcp
Run:  a2ensite whatsUp
Run:  service apache2 restart
