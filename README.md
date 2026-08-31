# whatsUp
Comprehensive tool for scheduling nighttime observations

Installation:
1. cd /var/www/html
2. Run: git clone https://github.com/sifank/whatsUp.git
3. Move whatsUp.conf to /etc/apache2/sites-available
4. chown -R www-data:www-data /var/www/whatsUp
5. add "Listen 5004" to /etc/apache2/ports.conf
6. sudo ufw allow 5004/tcp
7. Run:  a2ensite whatsUp
8. Run:  service apache2 restart

Modify whatsUp.py, modify/add:
1. TELESCOPE [name, focal length]
2. CAMERAS [name, pixel size, width px, height px)
   
