# post_spider.py
# 抓取一个帖子中的所有楼层内容
import requests
from pyquery import PyQuery as pq

class post_spider(object):
    def __init__(self):
        object.__init__(self)
        self.session = requests.session()
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.116 Safari/537.36'
        }
        self.session.headers.update(headers)
        self.fileName = ''

    def crawlPost(self, url):
        r = self.session.get(url)
        doc = pq(r.text)
        bar = str(doc('a.card_title_fname').html()).replace(' ', '')
        title = str(doc('h3.core_title_txt').html())
        # create a file and write bar name and post name
        self.fileName = bar + '_' + title + '.txt'
        try:
            file = open(self.fileName, 'w', encoding='utf-8')
            file.truncate()
            file.write(bar + '\n')
            file.write(title + '\n\n')
            # crawl the post and write every floor into the file
            pn = int(doc('li.l_reply_num span.red:last-child').html())
            for n in range(pn):
                urll = url + '?pn=' + str(n+1)
                print('crawling page:' + urll)
                self.crawlPage(urll, file)

            file.close()
        except FileNotFoundError as e:
            print('Except:', e)


    def crawlPage(self, url, file):
        floorList = []
        floorDict = {}
        r = self.session.get(url)
        # test crawing web html
        # print(r.text)
        # f = open('testweb.html', 'w', encoding='utf-8')
        # f.write(r.text)
        # f.close()

        doc = pq(r.text)
        # 将每一楼的html存到floorList中
        doc('.l_post.l_post_bright.j_l_post').each(lambda index:floorList.append(doc(this).html()))

        for item in floorList:
            doc = pq(item)
            floorDict['author'] = doc('.p_author_name.j_user_card').text()
            if(doc('.louzhubiaoshi_wrap')):
                floorDict['louzhu'] = 1
            else:
                floorDict['louzhu'] = 0
            floorDict['content'] = doc('.d_post_content.j_d_post_content').text()
            floorDict['time'] = doc('div.post-tail-wrap span:nth-child(4)').text()

            # write into file
            file.write(floorDict['author'] + '\n')
            if(floorDict['louzhu'] == 1):
                file.write('[LZ]\n')
            file.write(str(floorDict['content']) + '\n')
            file.write(floorDict['time'] + '\n')
            file.write('---------------------------------------------------------\n')

    def printAFloor(self, floorDict):
        file = open(self.fileName, 'w', encoding='utf-8')
        file.write(floorDict['author'] + '\n')
        if(floorDict['louzhu'] == 1):
            file.write('[LZ]\n')
        file.write(str(floorDict['content']) + '\n')
        file.write(floorDict['time'] + '\n')


if __name__ == '__main__':
    s = post_spider()
    s.crawlPost('http://tieba.baidu.com/p/4667218052')      # 念娇奴
    s.crawlPost('http://tieba.baidu.com/p/4643598682')      # 潇湘溪苑
    s.crawlPost('http://tieba.baidu.com/p/3765836395')      # 潇湘同城
    s.crawlPost('http://tieba.baidu.com/p/4643635451')      # 哥哥
    s.crawlPost('http://tieba.baidu.com/p/4642837650')      # 找妹妹
    s.crawlPost('http://tieba.baidu.com/p/4643915470')      # 找哥哥
    s.crawlPost('http://tieba.baidu.com/p/4184969803')      # 小贝
    s.crawlPost('http://tieba.baidu.com/p/4674929606')      # 惩戒
