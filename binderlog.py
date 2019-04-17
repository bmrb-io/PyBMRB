import requests
import json
import sys


def binderlog(date):
    c=0
    response = requests.get("https://archive.analytics.mybinder.org/events-{}.jsonl".format(date))
    #https://archive.analytics.mybinder.org/index.jsonl
    try:
        data = [json.loads(l) for l in response.iter_lines()]
        for i in data:
            if 'PyBMRB' in i['spec'].split("/"):
                #print (i)
                c+=1
    except json.decoder.JSONDecodeError:
        c = -1

    return c

def month_log(y,m):
    c=0
    for i in range(1,32):
        date = '{}-{:02d}-{:02d}'.format(y,m,i)
        n=binderlog(date)
        if n > 0:
            c+=n
            print ("{} : {}".format(date,n))
    print ("   {}-{:02d} : {}".format(y,m,c))

if __name__=="__main__":
    y=int(sys.argv[1])
    m=int(sys.argv[2])
    month_log(y,m)