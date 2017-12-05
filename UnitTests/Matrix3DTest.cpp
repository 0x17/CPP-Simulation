
#include "Matrix3DTest.h"
#include <array>

void assertEqualsExampleMatrix(Matrix3D<int> &m) {
	for(int h=0; h < 1; h++) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				ASSERT_EQ(i * 4 + j + 1, m.at(h, i, j));
			}
		}
	}
}

TEST_F(Matrix3DTest, testNestedVectorConstructor) {
	assertEqualsExampleMatrix(*m);
}

TEST_F(Matrix3DTest, testLambdaConstructor) {
	Matrix3D<int> mx(1, 3, 4, [](int h, int i, int j) { return i * 4 + j + 1; });
	TestHelpers::matrixEquals(*m, mx);
}

TEST_F(Matrix3DTest, testCopyConstructor) {
	Matrix3D<int> mx(*m);
	assertEqualsExampleMatrix(mx);
}

TEST_F(Matrix3DTest, testDefaultConstructor) {
	Matrix3D<int> mx;
    ASSERT_EQ(0, mx.getL());
	ASSERT_EQ(0, mx.getM());
	ASSERT_EQ(0, mx.getN());
}

TEST_F(Matrix3DTest, testFixedSizeEmptyConstructor) {
	Matrix3D<int> mx(1, 2, 3);
    ASSERT_EQ(1, mx.getL());
	ASSERT_EQ(2, mx.getM());
	ASSERT_EQ(3, mx.getN());
}

TEST_F(Matrix3DTest, testSingleValueConstructor) {
	Matrix3D<int> mx(1, 2, 2, 23);
    for(int h=0; h < 1; h++)
	    for (int i = 0; i < 2; i++)
		    for (int j = 0; j < 2; j++)
			    ASSERT_EQ(23, mx(h, i, j));
}

TEST_F(Matrix3DTest, testGetMN) {
    ASSERT_EQ(1, m->getL());
	ASSERT_EQ(3, m->getM());
	ASSERT_EQ(4, m->getN());
}

TEST_F(Matrix3DTest, testSubscriptOperator) {
    for (int h = 0; h<1; h++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                ASSERT_EQ(i * 4 + j + 1, (*m)(h, i, j));
            }
        }
    }
}

TEST_F(Matrix3DTest, testAssignOperator) {
	Matrix3D<int> mx = *m;
	assertEqualsExampleMatrix(mx);
}

TEST_F(Matrix3DTest, testResize) {
	m->resize(1, 1, 1);
    ASSERT_EQ(1, m->getL());
	ASSERT_EQ(1, m->getM());
	ASSERT_EQ(1, m->getN());
	ASSERT_EQ(1, m->at(0, 0, 0));
}

TEST_F(Matrix3DTest, testRow) {
	std::vector<int> firstRow = { 1, 2, 3, 4 };
	std::vector<int> secondRow = { 5, 6, 7, 8 };
	std::vector<int> thirdRow = { 9, 10, 11, 12 };
	std::array<std::vector<int>, 3> expRows = { firstRow, secondRow, thirdRow };

	std::vector<int> actualRow;
	for(int i=0; i<3; i++) {
		actualRow = m->row(0, i);
		TestHelpers::arrayEquals(expRows[i], actualRow);
	}
}

TEST_F(Matrix3DTest, testColumn) {
	std::vector<int> firstCol = { 1, 5, 9 };
	std::vector<int> secondCol = { 2, 6, 10};
	std::vector<int> thirdCol = { 3, 7, 11 };
	std::vector<int> fourthCol = { 4, 8, 12 };
	std::array<std::vector<int>, 4> expColumns = { firstCol, secondCol, thirdCol, fourthCol };

	std::vector<int> actualColumn;
	for (int i = 0; i<4; i++) {
		actualColumn = m->column(0, i);
		TestHelpers::arrayEquals(expColumns[i], actualColumn);
	}
}

TEST_F(Matrix3DTest, testForeach) {
	Matrix3D<char> entryVisited(1, 3, 4, 0);
	m->foreach([&entryVisited](int h, int i, int j, int mhij) { ASSERT_EQ(i * 4 + j + 1, mhij); entryVisited(h, i, j) = 1; });
    for (int h = 0; h<1; h++)
	    for (int i = 0; i<3; i++)
		    for (int j = 0; j<4; j++)
			    ASSERT_EQ(1, entryVisited(h, i, j));
}

TEST_F(Matrix3DTest, testForeach2) {
	Matrix3D<char> entryVisited(1, 3, 4, 0);
	m->foreach2([&entryVisited](int h, int i, int j) {  entryVisited(h, i, j) = 1; });
    for (int h = 0; h<1; h++)
	    for(int i=0; i<3; i++)
		    for(int j=0; j<4; j++)
			    ASSERT_EQ(1, entryVisited(h, i, j));
}

TEST_F(Matrix3DTest, testForeachAssign) {
	Matrix3D<char> entryVisited(1, 3, 4, 0);
	entryVisited.foreachAssign([](int h, int i, int j) {  return 1; });
    for(int h=0; h<1; h++)
	    for (int i = 0; i<3; i++)
		    for (int j = 0; j<4; j++)
			    ASSERT_EQ(1, entryVisited(h, i, j));
}

TEST_F(Matrix3DTest, testToString) {
	ASSERT_EQ("Matrix3D(l=1,m=3,n=4,\n{{{1,2,3,4},\n{5,6,7,8},\n{9,10,11,12}\n}\n}\n", m->toString());
}