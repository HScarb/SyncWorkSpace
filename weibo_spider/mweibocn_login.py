import requests
import re
import json
import base64
import time
import math
import random
from PIL import Image
from urllib.parse import quote_plus

# 构造 Request headers
agent = 'Mozilla/5.0 (Windows NT 6.2; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/49.0.2623.110 Safari/537.36'
global headers
headers = {
    "Host": "passport.weibo.cn",
    "Connection": "keep-alive",
    "Upgrade-Insecure-Requests": "1",
    'User-Agent': agent
}

session = requests.session()
# 访问登陆的初始页面
index_url = "https://passport.weibo.cn/signin/login"
session.get(index_url, headers=headers)

def get_su(username):
    '''
    对 email 地址和手机号码 先在javascript中encodeURIComponent
    对应Python3 中的是url.parse.quote_plus
    然后在base64加密后decode
    :param username:
    :return:
    '''
    username_quote = quote_plus(username)
    username_base64 = base64.b64encode(username_quote.encode('utf-8'))
    print(username_quote)
    print(username_base64.decode('utf-8'))
    return username_base64.decode('utf-8')

def login_pre(username):
    '''
    采用构造参数的方式
    :param username:
    :return:
    '''
    params = {
        'checkpin': '1',
        'entry': 'mweibo',
        'su': get_su(username),
        'callback': 'jsoncallback' + str(int(time.time() * 1000) + math.floor(random.random() * 100000))
    }

    pre_url = 'https://login.sina.com.cn/sso/prelogin.php'  # return a json {"retcode":-2001,"msg":"entry is empty","exectime":1}
    headers['Host'] = 'login.sina.com.cn'
    headers['Referer'] = index_url
    pre = session.get(pre_url, params = params, headers = headers)
    pa = r'\((.*?)\)'
    res = re.findall(pa, pre.text)
    print(res)
    if res == []:
        print('some error occurs, check your internet connection or your account.')
    else:
        js = json.loads(res[0])
        if js['showpin'] == 1:
            headers['Host'] = 'passport.weibo.cn'
            capt = session.get('https://passport.weibo.cn/captcha/image', headers = headers)
            capt_json = capt.json()
            capt_base64 = capt_json['data']['image'].split('base64')[1]
            with open('capt.jpg', 'wb') as f:
                f.write(base64.b64encode(capt_base64))
                f.close()
            im = Image.open('capt.jpg')
            im.show()
            im.close()
            cha_code = input('input captcha:\n>')
            return cha_code, capt_json['data']['pcid']
        else:
            return ''

def login(username, password, pincode):
    postdata = {
        'username': username,
        'password': password,
        'savestate': '1',
        'ec': '0',
        'pagerefer': '',
        "entry": "mweibo",
        "wentry": "",
        "loginfrom": "",
        "client_id": "",
        "code": "",
        "qq": "",
        "hff": "",
        "hfp": "",
    }
    if pincode == '':
        pass
    else:
        postdata['pincode'] = pincode[0]
        postdata['pcid'] = pincode[1]
    headers['Host'] = 'passport.weibo.cn'
    headers['Reference'] = index_url
    headers['Origin'] = 'https://passport.weibo.cn'
    headers['Content-Type'] = 'application/x-www-form-urlencoded'

    post_url = 'https://passport.weibo.cn/sso/login'
    login = session.post(post_url, data=postdata, headers=headers)
    print(login.cookies)
    print(login.status_code)
    js = login.json()
    print(js)
    uid = js['data']['uid']
    crossdomain = js['data']['crossdomainlist']
    cn = 'https:' + crossdomain['sina.com.cn']
    #
    headers['host'] = 'login.sina.com.cn'
    session.get(cn, headers=headers)
    headers['host'] = 'weibo.cn'
    ht = session.get('http://weibo.cn/%s/info' % uid, headers=headers)
    #
    print(ht.text)
    pa = r'<title>(.*?)</title>'
    res = re.findall(pa, ht.text)
    print(res)
    print('Hello %s !' %res[0])


if __name__ == '__main__':
    username = '13336008558'
    password = '621374si0'
    pincode = login_pre(username)
    login(username, password, pincode)
