// JewelGrid.h
#pragma once
#include "cocos2d.h"
#include "Jewel.h"

class JewelGrid : public cocos2d::Node
{
public:
	bool init(int col, int row);
	static JewelGrid * create(int col, int row);

	void changeMap();
	bool isDeadMap();

private:
	Jewel* createAJewel(int col, int row);
	bool isLegalJewel(Jewel * jewel, int x, int y);
	void setJewelPixPos(Jewel * jewel, float x, float y);

	void swapJewels(Jewel * jewel1, Jewel * jewel2);
	void swapJewelMove(Jewel * jewel);
	
private:
	virtual bool onTouchBegan(cocos2d::Touch * unused_touch, cocos2d::Event * unused_event);
	virtual void onTouchMoved(cocos2d::Touch * unused_touch, cocos2d::Event * unused_event);

private:
	// game logic
	bool canCrush();
	void doCrush();

	void onJewelsSwap(float dt);
	void onJewelsSwapBack(float dt);
	void onJewelRefresh(float dt);
	void onJewelCrush(float dt);

private:
	CC_SYNTHESIZE(int, m_row, Row);
	CC_SYNTHESIZE(int, m_col, Col);

	Jewel * m_jewelSelected;
	Jewel * m_jewelToSwap;
	
	std::vector<std::vector<Jewel*>> m_JewelBox;
	cocos2d::Vector<Jewel*> m_crushJewelBox;
	cocos2d::Vector<Jewel*> m_newJewelBox;
};