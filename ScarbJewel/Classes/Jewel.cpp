// Jewel.cpp
#include "Jewel.h"
USING_NS_CC;

bool Jewel::initByType(int x, int y, int type)
{
	if (!Sprite::init())
		return false;
	
	m_x = x;
	m_y = y;
	m_type = type;
	isSwapping = false;
	isCrushing = false;

	// set texture
	char name[100] = { 0 };
	sprintf(name, "jewel%d.png", m_type);
	this->initWithTexture(TextureCache::getInstance()->getTextureForKey(name));
	
	this->setAnchorPoint(Vec2(0, 0));

	return true;
}

Jewel* Jewel::createByType(int x, int y, int type)
{
	auto jewel = new Jewel();
	if(jewel && jewel->init())
	{
		jewel->autorelease();
		return jewel;
	}
	else
	{
		CC_SAFE_DELETE(jewel);
		return nullptr;
	}
}

void Jewel::crush()
{
	isCrushing = true;
	auto action = FadeOut::create(0.2);
	auto call = CallFunc::create([this]()
	{
		this->removeFromParent();
		isCrushing = false;
	});
	this->runAction(Sequence::create(action, call, nullptr));
}