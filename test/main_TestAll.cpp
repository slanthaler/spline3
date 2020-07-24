#include "gtest/gtest.h"

#include "spline3_test.hpp"
#include "tridiagonal_matrix_test.hpp"

int main(int argc, char **argv){
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
