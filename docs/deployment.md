# Deployment
This guide will instruct you on how to deploy this application on a apache2 server on a ubuntu machine. This guide comes with the requirement that you have the main application running with the docker compose solution
## Installation of software
First install the latest version of apache2:
```bash
sudo apt update
sudo apt install apache2
```

Make sure the system process of your newly installed apache server is running by checking with: 
```bash
sudo systemctl status apache2
```
If it's not running, run following command:
```bash
sudo systemctl start apache2
```

Then when it's running check your machines' IP adress with:
```bash
hostname -I
```
and navigate to said IP adress like http://<ip-adress> to see if you can access it.

## Configuration of apache

You will need to install two packages: proxy and proxy_http:
```bash
sudo a2enmod proxy
sudo a2enmode proxy_http
```
Then create a file like following:
```bash
sudo echo """<VirtualHost *:80>
    ServerName example.com
    ProxyPass /samplesheet http://localhost:8082/samplesheet
    ProxyPassReverse /samplesheet http://localhost:8082/samplesheet
</VirtualHost>
""" > /etc/apache2/sites-available/samplesheet.conf
```
Then you should be able to navigate to http://<ip-adress>/samplesheet and find the website!