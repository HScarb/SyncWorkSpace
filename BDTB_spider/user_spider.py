import requests
from pyquery import PyQuery as pq
import re
import urllib.parse

def print_list(list):   # 用于测试，将list中的数据分行打印
    for item in list:
        print(item)

class user_spider(object):
    def __init__(self):
        object.__init__(self)
        self.session = requests.Session()
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.116 Safari/537.36'
        }

    def crawl(self, ID):
        url = 'https://www.baidu.com/p/%s?from=tieba' % urllib.parse.quote(ID)
        detailurl = 'https://www.baidu.com/p/%s/detail' % urllib.parse.quote(ID)
        print(url)
        print(detailurl)

        r = requests.get(url)
        print(r.encoding)
        t = r.text
        # 将乱码转换成正确的字符串
        t = t.encode('latin1').decode('unicode_escape')
        t = t.encode('latin1').decode('utf-8')
        t = t.encode('utf-8').decode('utf-8')
        re_bars = re.compile(r'<ul class="honor-list clearfix .+<div>')
        print(t)
        # re_bar = re.compile(r'<li>.+<div>')
        result = re.findall(re_bars, t)
        print(result)
        doc = pq(result[0])
        print(doc.text())
        barList = re.findall(r'\d+?级 .+? ', doc.text())
        print_list(barList)

        r = requests.get(detailurl)
        t = r.text.encode('latin1').decode('utf-8')
        doc = pq(t)
        name = doc('a.link').text()
        print(name)
        profileList = []
        txt = doc('div.mod-profile').find('dl').text()
        print(txt)

        f = open(name + '.txt', 'w')
        f.truncate()
        f.write(name + '\n')
        f.write(url + '\n')
        f.write(detailurl + '\n')
        f.write(txt + '\n')
        for item in barList:
            f.write(item + '\n')
        f.close()


    def urllibcrawl(self):
        print()

if __name__ == '__main__':
    s = user_spider()
    s.crawl('李小呆萌萌哒')
