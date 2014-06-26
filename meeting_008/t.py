import schedule
import time
import os
from ftplib import FTP

def job():
    print("I'm working...")
    os.system('runipy -o plot_nrt_seaice.ipynb')
    os.system('ipython nbconvert --to html --template full_noinput.tpl plot_nrt_seaice.ipynb')

    ftp = FTP('ftp.zmaw.de')
    ftp.login(user='youruser', passwd='yourpass')
    ftp.cwd('outgoing/koldunov')
    ffile = open('plot_nrt_seaice.html','rb') 
    ftp.storbinary('STOR plot_nrt_seaice.html', ffile)
    ffile.close()
    ftp.quit()

schedule.every(30).minutes.do(job)
schedule.run_all()

while True:
    schedule.run_pending()
    time.sleep(1)

