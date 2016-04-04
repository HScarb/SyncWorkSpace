// JewelGrid.cpp
#include "JewelGrid.h"
USING_NS_CC;
const int FIRST_JEWEL_ID = 1;
const int LAST_JEWEL_ID = 6;

bool JewelGrid::init(int col, int row)
{
	if (!Node::init())
		return false;
	
	m_row = row;
	m_col = col;

	m_jewelSelected = nullptr;
	m_jewelToSwap = nullptr;

	m_JewelBox.resize(m_row);
	for(auto &vec : m_JewelBox)
	{
		vec.resize(m_col);
	}

	for (int x = 0; x < m_col;x++)
	{
		for (int y = 0; y < m_row;y++)
		{
			m_JewelBox[x][y] = createAJewel(x, y);
		}
	}

	while(isDeadMap())
	{
		changeMap();
	}

	// add touch listener
	auto listener = EventListenerTouchOneByOne::create();
	listener->setSwallowTouches(true);
	listener->onTouchBegan = CC_CALLBACK_2(JewelGrid::onTouchBegan, this);
	listener->onTouchMoved = CC_CALLBACK_2(JewelGrid::onTouchMoved, this);
	_eventDispatcher->addEventListenerWithSceneGraphPriority(listener, this);

	log("JewelGrid init!");
	return true;
}

JewelGrid* JewelGrid::create(int col, int row)
{
	JewelGrid * jewel = new JewelGrid();
	if(jewel && jewel->init(col, row))
	{
		jewel->autorelease();
	}
	else
	{
		CC_SAFE_DELETE(jewel);
		return nullptr;
	}
	return jewel;
}

void JewelGrid::changeMap()
{
	for (int x = 0; x < m_col;x++)
	{
		for (int y = 0; y < m_row;y++)
		{
			m_JewelBox[x][y]->removeFromParent();
			m_JewelBox[x][y] = createAJewel(x, y);
		}
	}
	log("update a new map!");
}

bool JewelGrid::isDeadMap()
{
	// simulate swap, and judge if it can be eliminated. If can't, this is a dead map
	auto swap = [](Jewel * jewel1, Jewel * jewel2)
	{
		auto temp = jewel1;
		jewel1 = jewel2;
		jewel2 = temp;
	};

	bool isDeadMap = true;

	for (int x = 0; x < m_col;x++)
	{
		for (int y = 0; y < m_row;y++)
		{
			// swap with left one
			if(x>0)
			{
				swap(m_JewelBox[x][y], m_JewelBox[x - 1][y]);
				if (canCrush())
					isDeadMap = false;
				swap(m_JewelBox[x][y], m_JewelBox[x - 1][y]);
			}
			// swap with right one
			if (x < m_col - 1)
			{
				swap(m_JewelBox[x][y], m_JewelBox[x + 1][y]);
				if (canCrush())
					isDeadMap = false;
				swap(m_JewelBox[x][y], m_JewelBox[x + 1][y]);
			}
			// swap with upper one
			if (y < m_row - 1)
			{
				swap(m_JewelBox[x][y], m_JewelBox[x][y + 1]);
				if (canCrush())
					isDeadMap = false;
				swap(m_JewelBox[x][y], m_JewelBox[x][y + 1]);
			}
			// swap with under one
			if(y<0)
			{
				swap(m_JewelBox[x][y], m_JewelBox[x][y - 1]);
				if (canCrush())
					isDeadMap = false;
				swap(m_JewelBox[x][y], m_JewelBox[x][y - 1]);
			}
			if (isDeadMap == false)
				break;
		}
	}
	// canCrush will store jewels in m_crushJewelBox, so clean it

	return isDeadMap;
}

Jewel* JewelGrid::createAJewel(int col, int row)
{
	Jewel * jewel = nullptr;

	while(1)
	{
		jewel = Jewel::createByType(col, row, random(FIRST_JEWEL_ID, LAST_JEWEL_ID));
		if (isLegalJewel(jewel))
			break;
	}
	setJewelPixPos(jewel, col, row);
}

bool JewelGrid::isLegalJewel(Jewel* jewel)
{
	bool isXLegal = true;
	bool isYLegal = true;

	if()
}

void JewelGrid::setJewelPixPos(Jewel* jewel, float x, float y)
{
}

void JewelGrid::swapJewels(Jewel* jewel1, Jewel* jewel2)
{
}

bool JewelGrid::onTouchBegan(cocos2d::Touch* unused_touch, cocos2d::Event* unused_event)
{
}

void JewelGrid::onTouchMoved(cocos2d::Touch* unused_touch, cocos2d::Event* unused_event)
{
}

bool JewelGrid::canCrush()
{
}

void JewelGrid::doCrush()
{
}

void JewelGrid::onJewelsSwap(float dt)
{
}

void JewelGrid::onJewelsSwapBack(float dt)
{
}

void JewelGrid::onJewelRefresh(float dt)
{
}

void JewelGrid::onJewelCrush(float dt)
{
}