import requests
from pyquery import PyQuery as pq
from bs4 import BeautifulSoup
import urllib
import urllib.parse
from user_spider import user_spider
from post_spider import post_spider

def print_list(list):   # 用于测试，将list中的数据分行打印
    for item in list:
        print(item)

class bar_spider(object):
    def __init__(self):
        object.__init__(self)
        self.session = requests.session()
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.116 Safari/537.36'
        }

    def crawlPage(self, url):
        r = requests.get(url)
        soup = BeautifulSoup(r.text, 'html.parser')
        postList = soup.find_all('a', class_ = 'j_th_tit ')
        # print_list(postList)
        ps = post_spider()
        for item in postList:
            ps.crawlPost('http://tieba.baidu.com' + item.get('href'))

    def crawl(self, barName, page):
        for x in range(page):
            pn = x * 50
            url = 'http://tieba.baidu.com/f?kw=' + urllib.parse.quote(barName) + '&pn=' + str(pn)
            print('Crawling page ' + str(x) + '...')
            self.crawlPage(url)


if __name__ == '__main__':
    s = bar_spider()
    s.crawl('脱团', 2)