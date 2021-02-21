import time
import urllib.parse
import urllib.request
from urllib.parse import urlparse
from bs4 import BeautifulSoup
from . import utils
import random
import re

def repeatmask(sequence,dnasource='mouse'):
    oksource = ['vertebrate','mammal','human','rodent','mouse','rat','danio','drosophila','elegans']
    assert dnasource in oksource
    params = {'sequence':sequence,
                'engine':'wublast',
                'speed':'default',
                'dnasource':dnasource,
                'ReturnFormat':'links',
                'ReturnMethod':'html'}

    baseURL = 'http://www.repeatmasker.org'
    actionURL = '/cgi-bin/WEBRepeatMasker'

    data = urllib.parse.urlencode(params)
    data = data.encode('ascii')

    req = urllib.request.Request(baseURL+actionURL, data)

    with urllib.request.urlopen(req) as response:
       html_text = response.read()

    soup = BeautifulSoup(html_text,'html.parser')
    resultLink = soup.find_all('a')[1].get('href')
    resultParse = urlparse(resultLink)
    resultLink = baseURL + resultParse.path
    utils.eprint(f'Results will be available at URL: {resultLink}\nParsing Now...')

    time.sleep(5)
    sleepInterval = 5
    maskDone = False
    termText = ['Masked File','No repetitive sequences were detected']
    while not maskDone:
        with urllib.request.urlopen(resultLink) as response:
            #print(response.getcode())
            result_text = response.read()
            #print(result_text)
        for term in termText:
            if term in str(result_text):
                maskDone = True
        utils.eprint(f'Results not ready. Trying again in {sleepInterval} seconds.')
        time.sleep(sleepInterval)

    # If no repeats, return orig sequence
    if termText[1] in str(result_text):
        return(sequence)
    else:
        res_soup = BeautifulSoup(result_text,'html.parser')
        summary_text = res_soup.find_all("pre")[1].text
        utils.eprint('*** SUMMARY OF REPEATMASKER JOB ***')
        utils.eprint(summary_text)

        files = []
        for link in res_soup.find_all('a'):
            files.append(link.get('href'))

        for file in files:
            #print(file)
            if 'masked' in file:
                #print("Got it!")
                with urllib.request.urlopen(baseURL+file) as response:
                    res = response.read()
                res = res.decode("utf-8")
                res = res[res.find('\n'):]
                res = res.replace("\n","")
                utils.eprint("Returning Masked Sequence")
                return(res)

def test():
    test_str_len = 30
    test_seq = ''.join(random.choices(['a','g','c','t'], k=test_str_len)) + 'agagagagagagagagagagagagagagagagagaga' + ''.join(random.choices(['a','g','c','t'], k=test_str_len))
    test_seq_random = ''.join(random.choices(['a','g','c','t'], k=int(test_str_len*2)))
    test_results = repeatmask(test_seq,dnasource='mouse')
    print(test_results)

if __name__ == "__main__":
    test()
