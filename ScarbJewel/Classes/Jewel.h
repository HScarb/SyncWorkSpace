// Jewel.h
#pragma once
#include "cocos2d.h"

class Jewel: public cocos2d::Sprite
{
public:
	bool initByType(int x, int y, int type);
	static Jewel * createByType(int x, int y, int type);

	void crush();

	CC_SYNTHESIZE(int, m_x, x);
	CC_SYNTHESIZE(int, m_y, y);
	CC_SYNTHESIZE(int, m_type, Type);
	CC_SYNTHESIZE(bool, isSwapping, IsSwapping);
	CC_SYNTHESIZE(bool, isCrushing, IsCrushing);
};