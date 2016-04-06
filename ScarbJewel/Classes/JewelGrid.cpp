// JewelGrid.cpp
#include "JewelGrid.h"
USING_NS_CC;
const int FIRST_JEWEL_ID = 1;
const int LAST_JEWEL_ID = 6;
const int GRID_WIDTH = 48;
const int MOVE_DURATION = 0.2;

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
		if (isLegalJewel(jewel, col, row))
			break;
	}
	setJewelPixPos(jewel, col, row);
}

bool JewelGrid::isLegalJewel(Jewel* jewel, int x, int y)
{
	bool isXLegal = true;
	bool isYLegal = true;

	if(jewel->getx() > 1)
	{
		if(jewel->getType() == m_JewelBox[x-1][y]->getType() && m_JewelBox[x-2][y]->getType())
		{
			isXLegal = false;
		}
	}
	if(y>1)
	{
		if(jewel->getType() == m_JewelBox[x][y-1] && m_JewelBox[x][y-2])
		{
			isYLegal = false;
		}
	}
	return isXLegal && isYLegal;
}

void JewelGrid::setJewelPixPos(Jewel* jewel, float x, float y)
{
	jewel->setPosition(x * GRID_WIDTH, y * GRID_WIDTH);
}

void JewelGrid::swapJewels(Jewel* jewel1, Jewel* jewel2)
{
	_eventDispatcher->pauseEventListenersForTarget(this);
	auto temp = m_JewelBox[jewel1->getx()][jewel1->gety()];
	m_JewelBox[jewel1->getx()][jewel1->gety()] = m_JewelBox[jewel2->getx()][jewel2->gety()];
	m_JewelBox[jewel2->getx()][jewel2->gety()] = temp;

	auto tempX = jewel1->getx();
	jewel1->setx(jewel2->getx());
	jewel2->setx(tempX);

	auto tempY = jewel1->gety();
	jewel1->setx(jewel2->getx());
	jewel2->setx(tempY);

	swapJewelMove(jewel1);
	swapJewelMove(jewel2);
}

// run move action for jewel, to move to a new position
void JewelGrid::swapJewelMove(Jewel* jewel)
{
	jewel->setIsSwapping(true);
	auto move = MoveTo::create(MOVE_DURATION, Vec2(jewel->getx()*GRID_WIDTH, jewel->gety() * GRID_WIDTH));
	auto call = CallFunc::create([jewel]()
	{
		jewel->setIsSwapping(false);
	});
	jewel->runAction(Sequence::create(move, call, nullptr));
}

bool JewelGrid::onTouchBegan(cocos2d::Touch* unused_touch, cocos2d::Event* unused_event)
{
	return true;
}

void JewelGrid::onTouchMoved(cocos2d::Touch* unused_touch, cocos2d::Event* unused_event)
{
}

// judge if the gird can be crushed
// if can, put jewels which can crush into crushedBox, wait to delete
bool JewelGrid::canCrush()
{
	int count = 0;
	Jewel * jewelBegin = nullptr;
	Jewel * jewelNext = nullptr;
	
	// traverse all jewels by column
	for (int x = 0; x < m_col;x++)
	{
		for (int y = 0; y < m_row;y++)
		{
			count = 1;
			jewelBegin = m_JewelBox[x][y];
			jewelNext = m_JewelBox[x][y + 1];

			while(jewelBegin->getType() == jewelNext->getType())
			{
				count++;
				int nextIndex = y + count;
				if (nextIndex > m_row - 1)
					break;
				jewelNext = m_JewelBox[x][nextIndex];
			}
			if(count >= 3)
			{
				for (int i = 0; i < count; i++)
				{
					auto jewel = m_JewelBox[x][y + i];
					m_crushJewelBox.pushBack(jewel);
				}
			}
			y += count;
		}
	}
	// traverse all jewels by row
	for (int y = 0; y < m_row;y++)
	{
		for (int x = 0; x < m_col;x++)
		{
			count = 1;
			jewelBegin = m_JewelBox[x][y];
			jewelNext = m_JewelBox[x + 1][y];

			while(jewelBegin->getType() == jewelNext->getType())
			{
				count++;
				int nextIndex = x + count;
				if (nextIndex > m_col - 1)
					break;
				jewelNext = m_JewelBox[nextIndex][y];
			}
			if(count >= 3)
			{
				for (int i = 0; i < count;i++)
				{
					auto jewel = m_JewelBox[x + i][y];
					m_crushJewelBox.pushBack(jewel);
				}
			}
			x += count;
		}
	}
	// if crushBox isn't empty, return true
	if (!m_crushJewelBox.empty())
		return true;
	else
		return false;
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