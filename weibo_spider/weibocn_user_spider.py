import requests
import json, re
import os
from bs4 import BeautifulSoup
from copy import deepcopy
from sys import stderr

class user_spider:
    def __init__(self, cookie,
                 max_num_weibo = 10,
                 max_num_fans = 10,
                 max_num_follow=10,
                 max_num_page=50,
                 wfilter='all',
                 return_type='string'):
        '''

        :param cookie: 登陆账户的cookie
        :param max_num_weibo: 最多抓取的微博条数
        :param max_num_fans: 最多抓取的分数数目
        :param max_num_follow: 最多抓取的关注
        :param max_num_page: 最多抓取的页面数
        :param wfilter: 抓取的微博类型 'all' 表示所有微博 'original'表示原创微博
        :param return_type: 返回值类型 'string'返回字符串化的json数据， 'json'返回一个json对象
        :return:
        '''
        self.INF = 10**9
        if max_num_weibo < 0:
            self.max_num_weibo = self.INF
        else:
            self.max_num_weibo = max_num_weibo
        if max_num_fans < 0:
            self.max_num_fans = self.INF
        else:
            self.max_num_fans = max_num_fans
        if max_num_follow < 0:
            self.max_num_follow = self.INF
        else:
            self.max_num_follow = max_num_follow
        if max_num_page < 0:
            self.max_num_page = self.INF
        else:
            self.max_num_page = max_num_page
        self.wfilter = wfilter
        self.data = None
        self.return_type = return_type

        # some constant, DO NOT CHANGE THEM
        RED_V, BLUE_V, VIP, ACTIVE = 'RED_V', 'BLUE_V', 'VIP', 'ACTIVE'
        self.verify_table = {'http://u1.sinaimg.cn/upload/2011/07/28/5338.gif': RED_V,
                        'http://u1.sinaimg.cn/upload/2011/07/28/5337.gif': BLUE_V,
                        'http://u1.sinaimg.cn/upload/h5/img/hyzs/donate_btn_s.png': VIP,
                        'http://u1.sinaimg.cn/upload/2011/08/16/5547.gif': ACTIVE}

        self.headers = {
                'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                'Accept-Encoding':'gzip, deflate, sdch',
                'Accept-Language':'zh-CN,zh;q=0.8,en;q=0.6,pl;q=0.4,zh-TW;q=0.2,ru;q=0.2',
                'Connection':'keep-alive',
                'Cookie':'',
                'Host':'weibo.cn',
                'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/43.0.2357.132 Safari/537.36'
            }
        self.ALL, self.ORIGINAL = 'all', 'original'
        # END CONSTANT DEFINITION

        if self.wfilter != self.ALL and self.wfilter != self.ORIGINAL:
            stderr.write('Invalid wfilter parameter, expect: all | original\n')
            exit(1)
        self.headers['Cookie'] = cookie

    def crawl(self, url='http://weibo.cn/yaochen'):
        self.data = {
            'url': '',
            'nickname': '',
            'verify_type': '',
            'verify_info': '',
            'num_weibo': -1,
            'num_fans': -1,
            'num_follow': -1,
            'gender': '',
            'birthday': '',
            'weibo': [],
            'fans': [],
            'follow': [],
            'location': '',
            'relationship_status': '',
            'sexual_orientation': '',
            'good_at': '',
            'self-intro': '',
            'tags': []
        }
        self.data['url'] = self.__normalize_url(url)

        # 1. crawl weibo text
        turl = self.data['url']
        if self.wfilter == self.ORIGINAL:
            turl += '?filter=1'
        req = self.__get_request(turl)
        soup = BeautifulSoup(req.text, 'lxml')
        if self.__abnormal(soup):
            stderr.write('Error: %s is abnormal user.\n' % url)
            return deepcopy(self.data)
        n_page = self.__parse_max_pages(soup)
        self.data['weibo'] = self.__parse_weibo(soup)
        self.data['verify_type'] = self.__parse_user_verify_type(soup)
        info_url = self.__parse_info_url(soup)
        fans_url, follow_url = self.__parse_fans_and_follow_url(soup)
        # page by page
        for i in range(2, 1 + n_page):
            if i > self.max_num_page or len(self.data['weibo']) >= self.max_num_weibo:
                break
            turl = self.data['url'] + '?page=' + str(i)
            if self.wfilter == self.ORIGINAL:
                turl = self.data['url'] + '?filter=1&page=' + str(i)
            req = self.__get_request(turl)
            soup = BeautifulSoup(req.text, 'lxml')
            self.data['weibo'] += self.__parse_weibo(soup)
        self.data['weibo'] = self.data['weibo'][:self.max_num_weibo]

        # 2. crawl user's personal information
        req = self.__get_request(info_url)
        soup = BeautifulSoup(req.text, 'lxml')
        self.__parse_info_list(soup)

        # 3. crawl user's fans list
        if self.max_num_fans > 0 and self.data['num_fans'] > 0:
            req = self.__get_request(fans_url)
            soup = BeautifulSoup(req.text, 'lxml')
            fans_max_page = self.__parse_max_pages(soup)
            self.data['fans'] += self.__parse_fans(soup)
            for i in range(2, 1 + fans_max_page):
                if i > self.max_num_page or len(self.data['fans']) >= self.max_num_fans:
                    break
                req = self.__get_request(fans_url + '?page=' + str(i))
                soup = BeautifulSoup(req.text, 'lxml')
                self.data['fans'] += self.__parse_fans(soup)
            self.data['fans'] = self.data['fans'][:self.max_num_fans]

        # 4. crawl user's follow list
        if self.max_num_follow > 0 and self.data['num_follow'] > 0:
            req = self.__get_request(follow_url)
            soup = BeautifulSoup(req.text, 'lxml')
            follow_max_page = self.__parse_max_pages(soup)
            self.data['follow'] += self.__parse_follow(soup)
            for i in range(2, 1 + follow_max_page):
                if i > self.max_num_page or len(self.data['follow']) >= self.max_num_follow:
                    break
                req = self.__get_request(follow_url + '?page' + str(i))
                soup = BeautifulSoup(req.text, 'lxml')
                self.data['follow'] += self.__parse_follow(soup)
            self.data['follow'] = self.data['follow'][:self.max_num_follow]

        # 5. return final result
            if self.return_type == 'string':
                return json.dumps(self.data, ensure_ascii=False, sort_keys=True, indent=4).encode('utf-8', 'replace')
            else:
                return deepcopy(self.data)

        

    def __normalize_url(self, url):
        url = self.__remove_qmark(url)
        url = url.strip('/')
        url = url.replace('weibo.com', 'weibo.cn')
        return url

    def __abnormal(self, soup):
        return soup.get_text().find('您当前访问的用户状态异常') >= 0

    def __parse_follow(self, soup):
        ret = []
        for table in soup.find_all('table'):
            cur = {'url': '', 'verify_type': '', 'nickname': '', 'num_fans': -1}
            cur['url'] = table.find_all('a')[0]['href']
            cur['nickname'] = table.find_all('a')[1].get_text()
            img = list(filter(lambda x: x['src'] in self.verify_table, table.find_all('img')))
            if len(img) > 0:
                cur['verify_type'] = self.verify_table[img[0]['src']]
            for e in table.find_all('td'):
                for ch in e.children:
                    try:
                        if ch.startwith('粉丝'):
                            cur['num_fans'] = int(ch[:ch.find('人')].strip('粉丝'))
                    except:
                        continue
            ret.append(cur)
            return ret

    def __remove_qmark(self, url):
        while url.find('?') >= 0:
            url = url[:url.find('?')]
        return url

    def __parse_max_pages(self, soup):
        '''
        return total number of pages of one's weibo
        :param soup:
        :return:
        '''
        try:
            res = int(soup.find(id='pagelist').get_text().split('/')[-1].strip('页'))
        except:
            return 1
        return res

    def __get_request(self, url, max_try = 3):
        stderr.write(url + '\n')
        cnt = 0
        while cnt < max_try:
            try:
                req = requests.get(url, headers=self.headers)
            except:
                cnt += 1
                continue
            if req.status_code != requests.codes.ok:
                break
            return req
        # should not reach here if everything is right
        stderr.write(json.dumps(self.data, ensure_ascii=False, sort_keys=True, indent=4).encode('utf-8', 'replace').decode('utf-8'))
        stderr.write('Error: %s\n' + url)
        exit(1)

    def __parse_weibo(self, soup):
        '''
        return all weibo in a page as a list
        :param soup:
        :return:
        '''
        ret = []
        for block in soup.find_all('div', 'c'):
            divs = block.find_all('div')
            if len(divs) == 1:
                text = self.__trim_like(divs[0].get_text())
                ret.append(text)
            elif len(divs) == 2 or len(divs) == 3:
                text = self.__trim_like(divs[0].get_text())
                comment = self.__trim_like(divs[-1].get_text())
                ret.append(text + comment)
            elif len(divs) == 0:
                continue
            else:
                stderr.write('Error: invalid weibo page')
                exit(1)
        return ret


    def __trim_like(self, s):
        pos = s.rfind('赞')      # 从后开始找
        if pos >= 0:
            return s[:pos]
        else:
            return s

    def __parse_user_verify_type(self, soup):
        r = soup.find_all('div', attrs={'class': 'ut'})[0].find_all('img', attrs={'alt': 'V'})
        if len(r) == 0:
            return ''
        src = r[0]['src']
        if src in self.verify_table:
            return self.verify_table[src]
        return ''

    def __parse_info_url(self, soup):
        table = soup.find_all('table')[0]
        for e in table.find_all('a'):
            if e.get_text() == '资料':
                assert(e['href'][0] == '/')
                ret = 'http://weibo.cn' + e['href']
                return self.__remove_qmark(ret)
        # should not reach here in normal case
        stderr.write('Error: can not find info tag.\n')
        exit(1)

    def __parse_fans_and_follow_url(self, soup):
        table = soup.find_all('div', attrs={'class': 'tip2'})[0]
        fans_url, follow_url = None, None
        for e in table.children:
            try:
                text = e.get_text()
                if text.startswith('微博'):
                    self.data['num_weibo'] = int(text.strip('微博|[|]'))
                elif text.startswith('关注'):
                    self.data['num_follow'] = int(text.strip('关注|[|]'))
                    follow_url = self.__remove_qmark('http://weibo.cn' + e['href'])
                elif text.startswith('粉丝'):
                    self.data['num_fans'] = int(text.strip('粉丝|[|]'))
                    fans_url = self.__remove_qmark('http://weibo.cn' + e['href'])
            except:
                continue
        print(fans_url)
        print(follow_url)
        assert(follow_url != None and fans_url != None)
        return fans_url, follow_url

    def __parse_info_list(self, soup):
        arr = soup.find_all('div')
        table = None
        for i in range(len(arr)):
            if arr[i]['class'][0] == 'tip' and arr[i].get_text() == '基本信息':
                assert(i + 1 < len(arr))
                table = arr[i + 1]
        assert(table != None)
        for c in table.children:
            try:
                if c['href'].find('keyword') > 0:
                    self.data['tags'].append(c.get_text())
            except:
                pass
            pos = c.find(':')
            if pos == None:
                continue
            if pos < 0 or pos + 1 == len(c):
                try:
                    pos = c.find('：')
                except:
                    continue
            if pos < 0 or pos + 1 == len(c):
                continue
            key, val = c[:pos], c[pos + 1:]
            assert(len(key) > 0 and len(val) > 0)
            if key == '昵称': self.data['nickname'] = val
            elif key == '性别': self.data['gender'] = val
            elif key == '地区': self.data['location'] = val
            elif key == '生日': self.data['birthday'] = val
            elif key == '简介': self.data['self-intro'] = val
            elif key == '性取向': self.data['sexual_orientation'] = val
            elif key == '认证' or key == u'认证信息': self.data['verify_info'] = val
            elif key == '感情状况': self.data['relationship_status'] = val
            elif key == '达人': self.data['good_at'] = val
            else: stderr.write('Not include attribute: %s %s\n' % (key, val))

    def __parse_fans(self, soup):
        res = []
        for follow in soup.find_all('table'):
            cur = {'url': '', 'verify_type': '', 'num_fans': -1, 'nickname': ''}
            person = follow.find_all('td')[1]
            img = list(filter(lambda x: x['src'] in self.verify_table, follow.find_all('img')))
            if len(img) > 0:
                cur['verify_type'] = self.verify_table[img[0]['src']]
            cur['url'] = person.find('a')['href']
            pos = person.get_text().rfind('粉丝')
            plain = person.get_text()
            cur['nickname'] = plain[0:pos]
            cur['num_fans'] = int(plain[pos + len('粉丝')].split('人')[0])
            res.append(cur)
        return res

if __name__ == '__main__':
    s = user_spider('_T_WM=06684d13fba19ff9d3e661b7c7b00450; ALF=1471853772; SCF=AgHFOygh8at-1tx8HOSeWDWg7jVZ45xBGdw_1-8lUmcvpojBRDZ29n0TY72Z5opEoDGt5y4WTG3S18zSGGtCuxg.; SUB=_2A256l0N1DeTxGedH6lQT8CrIyTiIHXVWeG09rDV6PUJbktBeLU-tkW1EWyYOzYm-Ltc6AQHZwIc4IMEJgw..; SUBP=0033WrSXqPxfM725Ws9jqgMF55529P9D9WFfKLEr-bqh24UjD8W-mq-K5JpX5o2p5NHD95Qp1K2ceo5XShzXWs4Dqcjs-cSyMJiT; SUHB=0WGtHT6L2fupgQ; SSOLoginState=1469264677; gsid_CTandWM=4u8xCpOz5otoB1tjhvZPi82uw6Q')
    print(s.crawl())