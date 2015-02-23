/***************************************************************************
 *            test_grid_set.cc
 *
 *
 *  Copyright  2008  Ivan S. Zapreev, Pieter Collins
 *            ivan.zapreev@gmail.com, pieter.collins@cwi.nl
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "config.h"

#include "bdd_set.h"
#include "logging.h"
#include "graphics.h"
#include "zonotope.h"
#include "formula.h"
#include "expression.h"
#include "function_set.h"
#include "set_checker.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

void test_constructors() {
    ARIADNE_PRINT_TEST_COMMENT("Testing constructors.");
    BDDTreeSet set1;
    ARIADNE_TEST_EQUAL(set1.dimension(),0);

    BDDTreeSet set2(3);
    ARIADNE_TEST_EQUAL(set2.grid(),Grid(3));
    ARIADNE_TEST_EQUAL(set2.root_cell_height(),0);
    ARIADNE_TEST_EQUAL(set2.root_cell_coordinates(),array<int>(3, 0,0,0));
    ARIADNE_TEST_EQUAL(set2.enabled_cells(),bddfalse);

    Grid grid(4, 1.25);
    BDDTreeSet set3(grid, true);
    ARIADNE_TEST_EQUAL(set3.grid(),grid);
    ARIADNE_TEST_EQUAL(set3.root_cell_height(),0);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(),array<int>(4, 0,0,0,0));
    ARIADNE_TEST_EQUAL(set3.enabled_cells(),bddtrue);
    
    // Test copy constructor
    set2 = set3;
    ARIADNE_TEST_EQUAL(set2,set3);

    // Test cloning operator
    BDDTreeSet* pset = set2.clone();
    ARIADNE_TEST_EQUAL(set2, *pset);

    // Test construction from a bdd
    bdd enabled_cells = bdd_ithvar(2);
    BDDTreeSet set4(grid, 1, array<int>(4, 1,0,0,0), enabled_cells);
    ARIADNE_TEST_EQUAL(set4.grid(),grid);
    ARIADNE_TEST_EQUAL(set4.root_cell_height(),1);
    ARIADNE_TEST_EQUAL(set4.root_cell_coordinates(),array<int>(4, 1,0,0,0));
    ARIADNE_TEST_EQUAL(set4.enabled_cells(),enabled_cells);
    
    // If the dimension of root_cell_coordinates and grid differs an exception must be thrown.
    ARIADNE_TEST_FAIL(set4 = BDDTreeSet(grid, 1, array<int>(3, 0,0,0), enabled_cells));
}

void test_properties_subdivisions() {
    ARIADNE_PRINT_TEST_COMMENT("Testing properties.");
    // Test empty
    BDDTreeSet set0;
    // The zero-dimensional set is always empty
    ARIADNE_TEST_ASSERT(set0.empty());
    BDDTreeSet set1(2, false);
    ARIADNE_TEST_ASSERT(set1.empty());
    bdd enabled_cells = bdd_nithvar(0) & (bdd_ithvar(2) | bdd_ithvar(3));
    BDDTreeSet set2(Grid(3), 4, array<int>(3, 1,0,0), enabled_cells);
    ARIADNE_TEST_ASSERT(!set2.empty());
    
    // Test size
    ARIADNE_TEST_EQUAL(set1.size(), 0);
    ARIADNE_TEST_EQUAL(set2.size(), 4);
    
    // Test mince and recombine
    set1.mince(3);
    ARIADNE_TEST_EQUAL(set1.size(), 0);
    set1.recombine();
    ARIADNE_TEST_EQUAL(set1.size(), 0);    
    BDDTreeSet set3(Grid(3), 0, array<int>(3, 1,0,0), enabled_cells);
    set3.mince(1);
    ARIADNE_TEST_EQUAL(set3.size(), 4);
    set3.mince(2);
    ARIADNE_TEST_EQUAL(set3.size(), 24);
    set3.mince(3);
    ARIADNE_TEST_EQUAL(set3.size(), 24*8);
    set3.recombine();
    ARIADNE_TEST_EQUAL(set3.size(), 4);
        
    // Measure is not implemented yet
    ARIADNE_TEST_FAIL(set2.measure());
    
    // Test root_cell
    // raise an error if the set is zero dimensional
    ARIADNE_TEST_FAIL(set0.root_cell());
    ARIADNE_TEST_EQUAL(set1.root_cell(), Box(2, 0.0,1.0, 0.0,1.0));
    enabled_cells = bdd_nithvar(0) & bdd_ithvar(1) & bdd_nithvar(2) & bdd_ithvar(3);    
    set0 = BDDTreeSet(Grid(2), 3, array<int>(2, 0,0), enabled_cells);
    ARIADNE_TEST_EQUAL(set0.root_cell(), Box(2, 0.0,2.0, -2.0,2.0));
    ARIADNE_TEST_EQUAL(set2.root_cell(), Box(3, 2.0,4.0, 0.0,2.0, -2.0,2.0));   
    
    // Test bounding box
    ARIADNE_TEST_ASSERT(set1.bounding_box().empty());
    ARIADNE_TEST_EQUAL(set0.bounding_box(), Box(2, 1.5,2.0, -2.0,-1.0));       
    ARIADNE_TEST_EQUAL(set2.bounding_box(), Box(3, 2.0,4.0, 0.0,2.0, -2.0,0.0));       
}

void test_predicates() {
    ARIADNE_PRINT_TEST_COMMENT("Testing predicates.");
    // Test subset and superset
    BDDTreeSet set0;
    BDDTreeSet set1(Grid(2), false);
    // testing w.r.t. a zero-dimensional set should raise an error
    ARIADNE_TEST_FAIL(subset(set0, set1));
    ARIADNE_TEST_FAIL(superset(set0, set1));
    // testing sets with different grids should raise an error
    BDDTreeSet set2(Grid(2, 1.25), true);
    ARIADNE_TEST_FAIL(subset(set1, set2));
    ARIADNE_TEST_FAIL(superset(set1, set2));
    
    // check test with an empty set
    bdd enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    ARIADNE_TEST_ASSERT(subset(set1,set2));
    ARIADNE_TEST_ASSERT(!superset(set1,set2));
    
    // Sets with different root cell coordinates are disjoint
    set1 = BDDTreeSet(Grid(2), true);
    ARIADNE_TEST_ASSERT(!subset(set1,set2));
    ARIADNE_TEST_ASSERT(!superset(set1,set2));
    
    // Sets with the same root coordinates
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), bdd_ithvar(0));
    ARIADNE_TEST_ASSERT(subset(set2,set1));
    ARIADNE_TEST_ASSERT(!superset(set2,set1));
    
    // Subset is a reflexive relation
    ARIADNE_TEST_ASSERT(subset(set2,set2));
    ARIADNE_TEST_ASSERT(superset(set2,set2));    
    
    // Test disjoint and overlap
    // testing w.r.t. a zero-dimensional set should raise an error
    ARIADNE_TEST_FAIL(disjoint(set0, set1));
    ARIADNE_TEST_FAIL(overlap(set0, set1));
     // testing sets with different grids should raise an error
    set2 = BDDTreeSet(Grid(2, 1.25), true);
    ARIADNE_TEST_FAIL(subset(set1, set2));
    ARIADNE_TEST_FAIL(superset(set1, set2));
    // check test with an empty set
    set1 = BDDTreeSet(Grid(2), false);
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    ARIADNE_TEST_ASSERT(disjoint(set1,set2));
    ARIADNE_TEST_ASSERT(!overlap(set1,set2));
    ARIADNE_TEST_ASSERT(!disjoint(set1,set1));
    ARIADNE_TEST_ASSERT(overlap(set1,set1));
    
    // Sets with different root cell coordinates are disjoint
    set1 = BDDTreeSet(Grid(2), true);
    ARIADNE_TEST_ASSERT(!subset(set1,set2));
    ARIADNE_TEST_ASSERT(!superset(set1,set2));

    // Sets with the same root coordinates
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), bdd_ithvar(0));
    ARIADNE_TEST_ASSERT(!disjoint(set1,set2));
    ARIADNE_TEST_ASSERT(overlap(set1,set2));
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), bdd_ithvar(0) & bdd_ithvar(1));
    ARIADNE_TEST_ASSERT(disjoint(set1,set2));
    ARIADNE_TEST_ASSERT(!overlap(set1,set2));
    
    // A set always overlaps itself
    ARIADNE_TEST_ASSERT(!disjoint(set2,set2));
    ARIADNE_TEST_ASSERT(overlap(set2,set2));    
    
    // Test subset inclusion with Box and ConstraintSet
    ConstraintSet cs0;
    RealVariable x("x");
    RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);
    RealExpression invf_x = 2.0*x - y;
    RealExpression invf_y = y;
	List<RealExpression> expr;
	expr.append(invf_x);
	expr.append(invf_y);
	VectorFunction invf(expr,varlist);
	Box dom(2, -1.0,1.0, -1.0,1.0);
    ConstraintSet cs1(invf, dom);    
    varlist.append(RealVariable("z"));
	invf = VectorFunction(expr,varlist);
    ConstraintSet cs2(invf, dom);        
    // create the full-space constraint set
    expr.clear();
    VectorFunction zerof(0,2);
    ConstraintSet cs3(zerof, Box(0));
    // testing w.r.t. a zero-dimensional set or box should raise an error
    ARIADNE_TEST_FAIL(set0.subset(Box(2, 0.0,1.0, 0.0,1.0)));
    ARIADNE_TEST_FAIL(set1.subset(Box(0)));
    ARIADNE_TEST_FAIL(covers(cs1, set0));    
     // testing w.rt. to a set with different dimension should raise an error
    ARIADNE_TEST_FAIL(set1.subset(Box(3)));
    ARIADNE_TEST_FAIL(covers(cs2, set1));
    // check test with an empty set or empty box
    set1 = BDDTreeSet(Grid(2), false);
    ARIADNE_TEST_ASSERT(definitely(set1.subset(Box(2, 0.0,1.0, 0.0,1.0))));
    Box ebx = Box::empty_box(2);
    ARIADNE_TEST_ASSERT(!possibly(set2.subset(ebx)));
    ARIADNE_TEST_ASSERT(definitely(covers(cs1, set1)));
    // test with the full-space constraint set
    ARIADNE_TEST_ASSERT(definitely(covers(cs3, set2)));
    // test with a general set
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) & bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    Box bx1(2, 2.25,5.0, 4.25,6.5);
    Box bx2(2, 2.25,5.0, 4.75,6.5);
    ARIADNE_TEST_ASSERT(definitely(set2.subset(bx1)));
    ARIADNE_TEST_ASSERT(!possibly(set2.subset(bx2)));
    enabled_cells = bdd_nithvar(0) & bdd_nithvar(1);;
    set1 = BDDTreeSet(Grid(2), 2, array<int>(2, 0,0), enabled_cells);    
    enabled_cells = enabled_cells & bdd_nithvar(2) & bdd_nithvar(3);
    set2 = BDDTreeSet(Grid(2), 2, array<int>(2, 0,0), enabled_cells);  
    // plot("test_bdd_set_set1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set1);
    enabled_cells = enabled_cells & bdd_nithvar(4);
    BDDTreeSet set3(Grid(2), 2, array<int>(2, 0,0), enabled_cells);    
    ARIADNE_TEST_ASSERT(!possibly(covers(cs1, set1)));
    ARIADNE_TEST_ASSERT(possibly(covers(cs1, set2)));
    ARIADNE_TEST_ASSERT(!definitely(covers(cs1, set2)));
    ARIADNE_TEST_ASSERT(definitely(covers(cs1, set3)));
    
    // Test superset with Box
    // testing w.r.t. a zero-dimensional set should raise an error
    ARIADNE_TEST_FAIL(set0.superset(Box(2, 0.0,1.0, 0.0,1.0)));
     // testing w.rt. to a box with different dimension should raise an error
    ARIADNE_TEST_FAIL(set1.superset(Box(3)));
    // check test with an empty set or empty box
    set1 = BDDTreeSet(Grid(2), false);
    ARIADNE_TEST_ASSERT(!possibly(set1.superset(Box(2, 0.0,1.0, 0.0,1.0))));
    ARIADNE_TEST_ASSERT(definitely(set2.superset(ebx)));
    // test with a general set
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    bx1 = Box(2, 2.25,5.0, 4.25,6.5);
    bx2 = Box(2, 2.6,2.9, 4.1,5.9);
    ARIADNE_TEST_ASSERT(!possibly(set2.superset(bx1)));
    ARIADNE_TEST_ASSERT(definitely(set2.superset(bx2)));

    // Test disjoint and overlap with Box and ConstraintSet
    // testing w.r.t. a zero-dimensional set or box should raise an error
    ARIADNE_TEST_FAIL(set0.disjoint(Box(2, 0.0,1.0, 0.0,1.0)));
    ARIADNE_TEST_FAIL(set1.disjoint(Box(0)));
    ARIADNE_TEST_FAIL(disjoint(cs1, set0));
    ARIADNE_TEST_FAIL(disjoint(cs0, set0));
    ARIADNE_TEST_FAIL(set0.overlaps(Box(2, 0.0,1.0, 0.0,1.0)));
    ARIADNE_TEST_FAIL(set1.overlaps(Box(0)));
    ARIADNE_TEST_FAIL(overlaps(cs1, set0));
    ARIADNE_TEST_FAIL(overlaps(cs0, set1));
     // testing w.rt. to a box with different dimension should raise an error
    ARIADNE_TEST_FAIL(set1.disjoint(Box(3)));
    ARIADNE_TEST_FAIL(set1.overlaps(Box(3)));
    ARIADNE_TEST_FAIL(disjoint(cs2, set1));
    ARIADNE_TEST_FAIL(overlaps(cs2, set1));
    // check test with an empty set or empty box
    set1 = BDDTreeSet(Grid(2), false);
    ARIADNE_TEST_ASSERT(definitely(set1.disjoint(Box(2, 0.0,1.0, 0.0,1.0))));
    ARIADNE_TEST_ASSERT(!possibly(set1.overlaps(Box(2, 0.0,1.0, 0.0,1.0))));
    ARIADNE_TEST_ASSERT(definitely(disjoint(cs1, set1)));
    ARIADNE_TEST_ASSERT(!possibly(overlaps(cs1, set1)));    
    ebx = Box::empty_box(2);
    ARIADNE_TEST_ASSERT(definitely(set2.disjoint(ebx)));
    ARIADNE_TEST_ASSERT(!possibly(set2.overlaps(ebx)));
    // test with the full-space constraint set
    ARIADNE_TEST_ASSERT(definitely(overlaps(cs3, set2)));
    ARIADNE_TEST_ASSERT(!possibly(disjoint(cs3, set2)));
    // test with a general set
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set2 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    bx1 = Box(2, 2.25,5.0, 4.25,6.5);
    bx2 = Box(2, 1.25,2.25, 1.75,4.25);
    ARIADNE_TEST_ASSERT(!possibly(set2.disjoint(bx1)));
    ARIADNE_TEST_ASSERT(definitely(set2.disjoint(bx2)));     
    ARIADNE_TEST_ASSERT(definitely(set2.overlaps(bx1)));
    ARIADNE_TEST_ASSERT(!possibly(set2.overlaps(bx2)));
    enabled_cells = bdd_ithvar(1) & bdd_nithvar(2) & bdd_nithvar(3);
    set1 = BDDTreeSet(Grid(2), 4, array<int>(2, 0,0), enabled_cells);    
    enabled_cells = enabled_cells & bdd_ithvar(4) & bdd_nithvar(5) & bdd_ithvar(6);
    set2 = BDDTreeSet(Grid(2), 4, array<int>(2, 0,0), enabled_cells);    
    enabled_cells = enabled_cells & bdd_nithvar(7);
    set3 = BDDTreeSet(Grid(2), 4, array<int>(2, 0,0), enabled_cells);    
    ARIADNE_TEST_ASSERT(definitely(overlaps(cs1, set1)));     
    ARIADNE_TEST_ASSERT(!possibly(disjoint(cs1, set1)));     
    ARIADNE_TEST_ASSERT(possibly(overlaps(cs1, set2)));     
    ARIADNE_TEST_ASSERT(possibly(disjoint(cs1, set2)));     
    ARIADNE_TEST_ASSERT(!possibly(overlaps(cs1, set3)));     
    ARIADNE_TEST_ASSERT(definitely(disjoint(cs1, set3)));     

}

void test_operations() {
    ARIADNE_PRINT_TEST_COMMENT("Testing operations.");
    
    // Test clear
    BDDTreeSet set0;
    // test clearing of a zero-dimensional set
    set0.clear();
    ARIADNE_TEST_ASSERT(set0.empty());
    // test clearing of a non-empty set
    BDDTreeSet set1(Grid(3, 1.25), 2, array<int>(3, 1,1,1), bdd_ithvar(1));
    ARIADNE_TEST_ASSERT(!set1.empty());
    set1.clear();
    ARIADNE_TEST_ASSERT(set1.empty());
    ARIADNE_TEST_EQUAL(set1.grid(), Grid(3, 1.25));
    ARIADNE_TEST_EQUAL(set1.root_cell_height(), 0);
    ARIADNE_TEST_EQUAL(set1.root_cell_coordinates(), array<int>(3, 0,0,0));
    ARIADNE_TEST_EQUAL(set1.enabled_cells(), bddfalse);
    
    // Test minimize_height
    // raise an error if the set is zero dimensional
    ARIADNE_TEST_FAIL(set0.minimize_height());
    set1 = BDDTreeSet(Grid(3), true);
    BDDTreeSet set2 = set1;
    // No changes if the height is zero.
    ARIADNE_TEST_CHECK(set1.minimize_height(), 0);
    ARIADNE_TEST_EQUAL(set1, set2);

    bdd enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set1 = BDDTreeSet(Grid(3), 4, array<int>(3, 1,1,1), enabled_cells);
    set2 = set1;
    ARIADNE_TEST_CHECK(set1.minimize_height(), 2);
    ARIADNE_TEST_ASSERT(set1 != set2);
    ARIADNE_TEST_EQUAL(set1.root_cell_height(), 2);
    ARIADNE_TEST_EQUAL(set1.root_cell_coordinates(), array<int>(3, 2,1,2));
    enabled_cells = bdd_ithvar(1) | bdd_ithvar(2);
    ARIADNE_TEST_EQUAL(set1.enabled_cells(), enabled_cells);
    ARIADNE_TEST_EQUAL(set1.root_cell(), Box(3, 2.0,3.0, 2.0,4.0, 4.0,6.0));
    
    // Test increase height
    // raise an error if the set is zero dimensional
    ARIADNE_TEST_FAIL(set0.increase_height(5));
    // increase is the opposite of minimze
    ARIADNE_TEST_CHECK(set1.increase_height(4), 4);
    ARIADNE_TEST_EQUAL(set1, set2);
    // if the new height is smaller than the current one the set is not changed
    ARIADNE_TEST_EQUAL(set1.increase_height(1), 4);
    ARIADNE_TEST_EQUAL(set1, set2);
    // test with arbitrary sets    
    BDDTreeSet set3 = BDDTreeSet(2, true);
    ARIADNE_TEST_CHECK(set3.increase_height(5), 5);
    ARIADNE_TEST_EQUAL(set3.grid(),Grid(2));
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), 5);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), array<int>(2, 0,0));
    ARIADNE_TEST_EQUAL(set3.root_cell(), Box(2, -2.0,2.0, -2.0,6.0));
    enabled_cells = bdd_nithvar(0) & bdd_ithvar(1) & bdd_ithvar(2) & bdd_nithvar(3) & bdd_nithvar(4);
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), enabled_cells);
    set3 = BDDTreeSet(Grid(2), 0, array<int>(2, -1,-1), bddtrue);
    ARIADNE_TEST_CHECK(set3.increase_height(6), 6);
    ARIADNE_TEST_EQUAL(set3.grid(),Grid(2));
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), 6);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), array<int>(2, 0,0));
    ARIADNE_TEST_EQUAL(set3.root_cell(), Box(2, -2.0,6.0, -2.0,6.0));
    enabled_cells = bdd_nithvar(0) & bdd_nithvar(1) & bdd_nithvar(2) & 
                    bdd_nithvar(3) & bdd_ithvar(4) & bdd_ithvar(5);
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), enabled_cells);
    
    // Test join
    // raise an error if one of the two sets is zero dimensional
    ARIADNE_TEST_FAIL(set0.adjoin(set1));
    ARIADNE_TEST_FAIL(join(set1, set0));
    // raise an error if the two sets have different dimension
    set2 = BDDTreeSet(Grid(2));
    ARIADNE_TEST_FAIL(set1.adjoin(set2));
    ARIADNE_TEST_FAIL(join(set1, set2));
    
    // joining an empty set do not change the set, except for minimization
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set1 = BDDTreeSet(Grid(3), 4, array<int>(3, 1,1,1), enabled_cells);
    set2 = BDDTreeSet(Grid(3), false);
    set3 = join(set1, set2);
    ARIADNE_TEST_ASSERT(set1 != set3);
    set1.minimize_height();
    ARIADNE_TEST_EQUAL(set1, set3);
    set2.adjoin(set1);
    ARIADNE_TEST_EQUAL(set1, set2);
    
    // test join of nonempty sets with the same root cell
    set1 = BDDTreeSet(Grid(3), 2, array<int>(3, 1,1,1), enabled_cells);
    bdd enabled_cells2 = bdd_nithvar(0) | (bdd_nithvar(1) & bdd_nithvar(3));
    set2 = BDDTreeSet(Grid(3), 2, array<int>(3, 1,1,1), enabled_cells2);
    set3 = join(set1, set2);
    ARIADNE_TEST_EQUAL(set3.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), set1.root_cell_height());
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), set1.root_cell_coordinates());
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), (enabled_cells | enabled_cells2));
    
    // test join of nonempty sets with different root cell
    enabled_cells = bdd_nithvar(0) | (bdd_nithvar(1) & bdd_nithvar(3));
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    set2 = BDDTreeSet(Grid(2), 0, array<int>(2, -1,-1), bddtrue);
    set3 = join(set1, set2);
    ARIADNE_TEST_EQUAL(set3.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), 6);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), array<int>(2, 0,0));
    enabled_cells = ( bdd_ithvar(0) & bdd_ithvar(1) & bdd_nithvar(2) &
                      (bdd_nithvar(3) | (bdd_nithvar(4) & bdd_nithvar(6))) )
                    |
                    ( bdd_nithvar(0) & bdd_nithvar(1) & bdd_nithvar(2) &
                      bdd_nithvar(3) & bdd_ithvar(4) & bdd_ithvar(5) ) ;
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), enabled_cells);
 
    // Test Intersection
    // raise an error if one of the two sets is zero dimensional
    ARIADNE_TEST_FAIL(set0.restrict(set1));
    ARIADNE_TEST_FAIL(intersection(set1, set0));
    // raise an error if the two sets have different dimension
    set2 = BDDTreeSet(Grid(4));
    ARIADNE_TEST_FAIL(set1.restrict(set2));
    ARIADNE_TEST_FAIL(intersection(set1, set2));
    
    // Intersection with an empty set make the set empty
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set1 = BDDTreeSet(Grid(3), 4, array<int>(3, 1,1,1), enabled_cells);
    set2 = BDDTreeSet(Grid(3), false);
    set3 = intersection(set1, set2);
    ARIADNE_TEST_ASSERT(set3.empty());

    // intersection with a superset set do not change the set, except for minimization
    set2 = BDDTreeSet(Grid(3), 4, array<int>(3, 1,1,1), bddtrue);
    set3 = intersection(set1, set2);
    ARIADNE_TEST_ASSERT(set1 != set3);
    set1.minimize_height();
    ARIADNE_TEST_EQUAL(set1, set3);
    set2.restrict(set1);
    ARIADNE_TEST_EQUAL(set1, set2);
    
    // test intersection of nonempty sets with the same root cell
    enabled_cells = bdd_ithvar(1) | (bdd_nithvar(1) & bdd_nithvar(3));
    set1 = BDDTreeSet(Grid(3), 2, array<int>(3, 1,1,1), enabled_cells);
    enabled_cells2 = bdd_nithvar(0) | (bdd_nithvar(1) & bdd_nithvar(3));
    set2 = BDDTreeSet(Grid(3), 2, array<int>(3, 1,1,1), enabled_cells2);
    set3 = intersection(set1, set2);
    ARIADNE_TEST_EQUAL(set3.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), set1.root_cell_height());
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), set1.root_cell_coordinates());
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), (enabled_cells & enabled_cells2));
    
    // intersection of nonempty sets with disjoint root cells must be empty
    enabled_cells = bdd_nithvar(0) | (bdd_nithvar(1) & bdd_nithvar(3));
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    set2 = BDDTreeSet(Grid(2), 0, array<int>(2, -1,-1), bddtrue);
    set3 = intersection(set1, set2);
    ARIADNE_TEST_ASSERT(set3.empty());
    
    // test intersection of partially overlapping sets
    set2 = BDDTreeSet(Grid(2), 1, array<int>(2, 2,2), bddtrue);
    set3 = intersection(set1, set2);
    ARIADNE_TEST_EQUAL(set3.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), 1);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), array<int>(2, 2,2));
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), bdd_nithvar(1));

    // Test Difference
    // raise an error if one of the two sets is zero dimensional
    ARIADNE_TEST_FAIL(set0.remove(set1));
    ARIADNE_TEST_FAIL(difference(set1, set0));
    // raise an error if the two sets have different dimension
    set2 = BDDTreeSet(Grid(4));
    ARIADNE_TEST_FAIL(set1.remove(set2));
    ARIADNE_TEST_FAIL(difference(set1, set2));

    // difference with an empty set do not change the set, except for minimization    
    enabled_cells = bdd_ithvar(0) & bdd_nithvar(1) & (bdd_ithvar(3) | bdd_ithvar(4));    
    set1 = BDDTreeSet(Grid(3), 4, array<int>(3, 1,1,1), enabled_cells);
    set2 = BDDTreeSet(Grid(3), false);
    set3 = difference(set1, set2);
    ARIADNE_TEST_ASSERT(set1 != set3);
    set1.minimize_height();
    ARIADNE_TEST_EQUAL(set1, set3);
    set1.remove(set2);
    ARIADNE_TEST_EQUAL(set1, set3);
    
    // test difference of nonempty sets with the same root cell
    enabled_cells = bdd_ithvar(0) | (bdd_nithvar(0) & bdd_nithvar(3));
    set1 = BDDTreeSet(Grid(3), 2, array<int>(3, 1,1,1), enabled_cells);
    enabled_cells2 = (bdd_nithvar(1) & bdd_nithvar(3));
    set2 = BDDTreeSet(Grid(3), 2, array<int>(3, 1,1,1), enabled_cells2);
    set3 = difference(set1, set2);
    ARIADNE_TEST_EQUAL(set3.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), set1.root_cell_height());
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), set1.root_cell_coordinates());
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), bdd_apply(enabled_cells, enabled_cells2, bddop_diff));
    
    // difference with a nonempty sets with disjoint root cell do no change the set
    enabled_cells = bdd_nithvar(0) | (bdd_nithvar(1) & bdd_nithvar(3));
    set1 = BDDTreeSet(Grid(2), 3, array<int>(2, 1,1), enabled_cells);
    set2 = BDDTreeSet(Grid(2), 0, array<int>(2, -1,-1), bddtrue);
    set3 = difference(set1, set2);
    ARIADNE_TEST_EQUAL(set3, set1);
    
    // test difference of partially overlapping sets
    set2 = BDDTreeSet(Grid(2), 1, array<int>(2, 2,2), bddtrue);
    set3 = difference(set1, set2);
    ARIADNE_TEST_EQUAL(set3.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set3.root_cell_height(), 2);
    ARIADNE_TEST_EQUAL(set3.root_cell_coordinates(), array<int>(2, 1,1));
    ARIADNE_TEST_EQUAL(set3.enabled_cells(), bddtrue);
 
}

void test_box_approximations() {
    ARIADNE_PRINT_TEST_COMMENT("Testing box approximations.");

    // Test adjoin over approximation of a box

    // raise an error if the set is zero-dimensional
    BDDTreeSet set0;
    ARIADNE_TEST_FAIL(set0.adjoin_over_approximation(Box(2, 1.0), 2));
    // raise an error if the dimensions are different
    BDDTreeSet set1(Grid(2), true);
    ARIADNE_TEST_FAIL(set1.adjoin_over_approximation(Box(3, 1.0), 2));
    // raise an error if the box is unbounded
    ARIADNE_TEST_FAIL(set1.adjoin_over_approximation(unbounded_box(2), 2));
    
    // adjoining an empty box do not change the set
    BDDTreeSet set2 = set1;
    set2.adjoin_over_approximation(Box::empty_box(2), 2);
    ARIADNE_TEST_EQUAL(set1, set2);
    
    // adjoining a proper subset of the current set do not change the set
    set2.adjoin_over_approximation(Box(2, 0.25,0.75, 0.25,0.75), 2);
    ARIADNE_TEST_EQUAL(set1, set2);
        
    // test adjoin with a general set
    bdd enabled_cells = bdd_nithvar(0) | (bdd_nithvar(1) & bdd_ithvar(3));
    set1 = BDDTreeSet(Grid(2, 2.0), 1, array<int>(2, 1,1), enabled_cells);
    set2 = set1;
    Box bx1 = Box(2, -1.9,2.9, 5.9,6.1);
    set2.adjoin_over_approximation(bx1, 1);

    // the original set is a subset of the new one
    ARIADNE_TEST_ASSERT(subset(set1, set2));
    // and not the other way around
    ARIADNE_TEST_ASSERT(!subset(set2, set1));
    // the box is a subset of the result
    ARIADNE_TEST_ASSERT(set2.superset(bx1));
    // but not a superset
    ARIADNE_TEST_ASSERT(!set2.subset(bx1));
    // compare with the correct result
    ARIADNE_TEST_EQUAL(set2.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set2.root_cell_height(), 4);
    ARIADNE_TEST_EQUAL(set2.root_cell_coordinates(), array<int>(2, 0,1));
    enabled_cells = ( bdd_nithvar(0) & bdd_nithvar(1) & bdd_ithvar(2) & 
                      ( (bdd_nithvar(3) & bdd_ithvar(5)) | (bdd_ithvar(3) & bdd_nithvar(5)) ) )
                    |
                    ( bdd_ithvar(0) & bdd_nithvar(1) & 
                      ( ( bdd_nithvar(2) & ((bdd_nithvar(3) & bdd_ithvar(5)) | (bdd_ithvar(3) & bdd_nithvar(5))) )
                        |
                        ( bdd_ithvar(2) & (bdd_nithvar(3) | (bdd_nithvar(4) & (bdd_nithvar(5) | bdd_ithvar(6)))) )
                    ) );
    ARIADNE_TEST_EQUAL(set2.enabled_cells(), enabled_cells);                          
   
    // adjoining with a greater depth gives a better result
    BDDTreeSet set3 = set1;
    set3.adjoin_over_approximation(bx1, 2);
    ARIADNE_TEST_ASSERT(subset(set1, set3));
    ARIADNE_TEST_ASSERT(subset(set3, set2));
    ARIADNE_TEST_ASSERT(!subset(set2, set3));
    set2 = set1;
    set3.adjoin_over_approximation(bx1, 3);
    ARIADNE_TEST_ASSERT(subset(set1, set2));
    ARIADNE_TEST_ASSERT(subset(set2, set3));
    ARIADNE_TEST_ASSERT(!subset(set3, set2));

    // Test adjoin with a set of zero measure
    bx1 = Box(2, 0.1,0.9, 1.0,1.0);
    set1 = BDDTreeSet(Grid(2));
    set2 = set1;
    set2.adjoin_over_approximation(bx1, 1);
    ARIADNE_TEST_ASSERT(set1 != set2);
    ARIADNE_TEST_EQUAL(set2.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set2.root_cell_height(), 1);
    ARIADNE_TEST_EQUAL(set2.root_cell_coordinates(), array<int>(2, 0,0));
    ARIADNE_TEST_EQUAL(set2.enabled_cells(), (bdd_ithvar(0) & bdd_nithvar(2)) | (bdd_nithvar(0) & bdd_ithvar(2)));                          

    set1 = BDDTreeSet(Grid(2, 2.0));
    set2 = set1;
    set2.adjoin_over_approximation(bx1, 1);
    ARIADNE_TEST_ASSERT(set1 != set2);
    ARIADNE_TEST_EQUAL(set2.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set2.root_cell_height(), 0);
    ARIADNE_TEST_EQUAL(set2.root_cell_coordinates(), array<int>(2, 0,0));
    ARIADNE_TEST_EQUAL(set2.enabled_cells(), bdd_nithvar(0));  
    
    bx1 = Box(2, 0.1,1.0, 0.5,0.5);
    set2 = set1;
    set2.adjoin_over_approximation(bx1, 1);
    ARIADNE_TEST_ASSERT(set1 != set2);
    ARIADNE_TEST_EQUAL(set2.grid(), set1.grid());
    ARIADNE_TEST_EQUAL(set2.root_cell_height(), 0);
    ARIADNE_TEST_EQUAL(set2.root_cell_coordinates(), array<int>(2, 0,0));
    ARIADNE_TEST_EQUAL(set2.enabled_cells(), bdd_nithvar(1));  


}

void test_set_approximations() {
    ARIADNE_PRINT_TEST_COMMENT("Testing set approximations.");

    // Test adjoin outer approximation of a set
/*
    Zonotope zt0(2);
    Zonotope zt1(2,2, 0.0,0.5,0.5,0.0, 0.0,0.0,1.0,0.0,0.0);
    
    // adjoin to an empty BDDTreeSet
    set1 = BDDTreeSet(Grid(2), false);
    set1.increase_height(zt1.bounding_box());
    ARIADNE_TEST_ASSERT(definitely(subset(zt1.bounding_box(), set1.root_cell())));
    Box dbox = set1.root_cell();
    plot("test_bdd_set_zt1",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),zt1);
    plot("test_bdd_set_zt1_bbox",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),zt1.bounding_box());
    
    ARIADNE_TEST_ASSERT(!definitely(zt1.disjoint(Box(2, 0.0,1.0, 0.0,2.0))));
    set1.adjoin_outer_approximation(zt1, 1);
    ARIADNE_TEST_ASSERT(!set1.empty());
    plot("test_bdd_set_zt1_adjoin",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set1);
*/

    RealVariable x("x");
    RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);
	RealExpression f_x = 0.5*x + 0.5*y;
	RealExpression f_y = y;
	List<RealExpression> expr;
	expr.append(f_x);
	expr.append(f_y);
	VectorFunction f(expr,varlist);
	Box dom(2, -1.0,1.0, -1.0,1.0);
	ImageSet is1(dom, f);
    RealExpression invf_x = 2.0*x - y;
    RealExpression invf_y = y;
    expr.clear();
	expr.append(invf_x);
	expr.append(invf_y);
	VectorFunction invf(expr,varlist);
    ConstraintSet cs1(invf, dom);
    // create the full-space constraint set
    expr.clear();
    VectorFunction zerof(0,2);
    ConstraintSet cs3(zerof, Box(0));
    
    // raise an error if the set is zero-dimensional
    BDDTreeSet set0;
    ARIADNE_TEST_FAIL(set0.adjoin_outer_approximation(is1, 1));
    ARIADNE_TEST_FAIL(set0.adjoin_lower_approximation(is1, 0, 1));
    ARIADNE_TEST_FAIL(set0.adjoin_inner_approximation(cs1, 0, 1));
    // raise an error if the dimensions are different
    BDDTreeSet set1o(Grid(4), true);
    ARIADNE_TEST_FAIL(set1o.adjoin_outer_approximation(is1, 1));
    ARIADNE_TEST_FAIL(set1o.adjoin_lower_approximation(is1, 1));
    ARIADNE_TEST_FAIL(set1o.adjoin_inner_approximation(cs1, Box(2, -1.0,1.0, -1.0,1.0), 1));

    // adjoin to an empty BDDTreeSet
    set1o = BDDTreeSet(Grid(2), false);
    set1o.increase_height(is1.bounding_box());
    ARIADNE_TEST_ASSERT(definitely(subset(is1.bounding_box(), set1o.root_cell())));
    Box dbox = set1o.root_cell();
    set1o.adjoin_outer_approximation(is1, 3);
    ARIADNE_TEST_ASSERT(!set1o.empty());
    plot("test_bdd_set_is1_outer",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set1o);
    BDDTreeSet set1l(Grid(2), false);
    set1l.adjoin_lower_approximation(is1, set1o.root_cell_height()/2 + 1, 3);
    ARIADNE_TEST_ASSERT(!set1l.empty());
    plot("test_bdd_set_is1_lower",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set1l);
    ARIADNE_TEST_ASSERT(definitely(subset(set1l, set1o)));
    ARIADNE_TEST_ASSERT(!possibly(superset(set1l, set1o)));
    BDDTreeSet set1i(Grid(2), false);
    set1i.adjoin_inner_approximation(cs1, set1o.root_cell_height()/2 + 1, 3);
    ARIADNE_TEST_ASSERT(!set1i.empty());
    plot("test_bdd_set_is1_inner",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set1i);
    ARIADNE_TEST_ASSERT(definitely(subset(set1i, set1l)));
    ARIADNE_TEST_ASSERT(!possibly(superset(set1i, set1l)));
    
    // adjoining an empty set do not change the set
    Box ebx = Box::empty_box(2);
    BDDTreeSet set2o = set1o;
    set2o.adjoin_outer_approximation(ebx, 5);
    ARIADNE_TEST_EQUAL(set1o, set2o);
    BDDTreeSet set2l = set1l;
    set2l.adjoin_lower_approximation(ebx, 5);
    ARIADNE_TEST_EQUAL(set1l, set2l);
    BDDTreeSet set2i = set1i;
    set2i.adjoin_inner_approximation(ebx, 0, 5);
    ARIADNE_TEST_EQUAL(set1i, set2i);

    // test adjoin with the full-space constraint set
    set2i = set1i;
    set2i.adjoin_inner_approximation(cs3, 0, 2);
    ARIADNE_TEST_EQUAL(set1i.root_cell(), set2i.root_cell());
    ARIADNE_TEST_EQUAL(set2i.enabled_cells(), bddtrue);
    
    // adjoin to a non-empty BDDTreeSet
	f_x = -0.25*x - 0.5*y;
	f_y = y;
	expr.clear();
	expr.append(f_x);
	expr.append(f_y);
	f = VectorFunction(expr,varlist);
	ImageSet is2(dom, f);
    invf_x = -4.0*x - 2.0*y;
    invf_y = y;
    expr.clear();
	expr.append(invf_x);
	expr.append(invf_y);
	invf = VectorFunction(expr,varlist);
    ConstraintSet cs2(invf, dom);
	
    set2o.adjoin_outer_approximation(is2, 4);
    plot("test_bdd_set_is1_is2_outer",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set2o);
    ARIADNE_TEST_ASSERT(definitely(subset(set1o, set2o)));
    ARIADNE_TEST_ASSERT(!possibly(subset(set2o, set1o)));
    set2l.adjoin_lower_approximation(is2, is2.bounding_box(), 4);
    plot("test_bdd_set_is1_is2_lower",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set2l);
    ARIADNE_TEST_ASSERT(definitely(subset(set1l, set2l)));
    ARIADNE_TEST_ASSERT(!possibly(subset(set2l, set1l)));
    set2i.adjoin_inner_approximation(cs2, is2.bounding_box(), 4);
    plot("test_bdd_set_is1_is2_inner",PlanarProjectionMap(2,0,1),dbox,Colour(1,0,1),set2i);
    ARIADNE_TEST_ASSERT(definitely(subset(set1i, set2i)));
    ARIADNE_TEST_ASSERT(!possibly(subset(set2i, set1i)));
    
    // Test adjoin with a point on a corner of the grid
    Box bx(2, 2.0,2.0, 0.0,0.0);
    set1o = BDDTreeSet(2);
    set1l = BDDTreeSet(2);
    set1i = BDDTreeSet(2);
    set1o.adjoin_outer_approximation(ImageSet(bx), 2);
    ARIADNE_TEST_EQUAL(set1o.size(), 4);
    set1l.adjoin_lower_approximation(bx, bx, 2);
    ARIADNE_TEST_ASSERT(set1l.empty());
    set1i.adjoin_inner_approximation(bx, bx, 2);
    ARIADNE_TEST_ASSERT(set1i.empty());
        	
}

void test_iterators_conversions_drawing() {
    ARIADNE_PRINT_TEST_COMMENT("Testing iterators.");
    
    // If the set is zero-dimensional an exception must be thrown 
    BDDTreeSet set0(0);
    BDDTreeSet::const_iterator it;
    ARIADNE_TEST_FAIL(it = set0.begin());

    // Test empty set iterator
    BDDTreeSet set1(2, false);
    for(it = set1.begin(); it != set1.end(); it++) {
        // the set is empty, never enter here
        ARIADNE_TEST_ASSERT(false);
    }

    // Test simple one-cell set
    set1 = BDDTreeSet(2, true);
    uint count = 0;
    for(it = set1.begin(); it != set1.end(); it++, count++) {
        // The set contains only one cell
        ARIADNE_TEST_EQUAL((*it), Box(2, 0.0,1.0, 0.0,1.0));
    }
    ARIADNE_TEST_EQUAL(count, 1);

    // Test complex set
    bdd enabled_cells = bdd_nithvar(0) & (bdd_ithvar(2) | bdd_ithvar(3));
    BDDTreeSet set2(Grid(2), 2, array<int>(2, 1,1), enabled_cells);
    array<Box> results(4);
    results[0] = Box(2, 2.0,2.5, 2.5,3.0);
    results[1] = Box(2, 2.5,3.0, 2.0,3.0);
    results[2] = Box(2, 2.0,2.5, 3.5,4.0);
    results[3] = Box(2, 2.5,3.0, 3.0,4.0);
    std::cout << "Start iteration." << std::endl;
    for(count = 0, it = set2.begin(); it != set2.end(); it++, count++) {
        std::cout << "Cell " << count << " = " << (*it) <<  std::endl;
        ARIADNE_TEST_EQUAL((*it), results[count]);
    }
    // The set contains 4 cells
    ARIADNE_TEST_EQUAL(count, 4);

    // Test conversion to a ListSet of Boxes
    ListSet<Box> boxlist(2);
    boxlist.push_back(results[0]);
    boxlist.push_back(results[1]);
    boxlist.push_back(results[2]);
    boxlist.push_back(results[3]);
    ListSet<Box> reslist = set2;
    ARIADNE_TEST_EQUAL(reslist, boxlist);
    
    // Test drawing
    plot("test_bdd_set_draw",PlanarProjectionMap(2,0,1),set2.bounding_box(),Colour(1,0,1),set2);

}

void test_projection() {
    ARIADNE_PRINT_TEST_COMMENT("Testing projections.");
    // Project down a zero-dimensional set should raise an error
    BDDTreeSet set1(0);
    BDDTreeSet set2;
    Vector<uint> indices(3);
    indices[0] = 1;
    indices[1] = 2;
    indices[2] = 4;
    ARIADNE_TEST_FAIL(set2 = project_down(set1, indices));
    // project with incorrect indices should raise an error
    set1 = BDDTreeSet(2);
    ARIADNE_TEST_FAIL(set2 = project_down(set1, indices));
    set1 = BDDTreeSet(4);
    ARIADNE_TEST_FAIL(set2 = project_down(set1, indices));
    // test projection of a real set.
    bdd enabled_cells = bdd_nithvar(0) & (bdd_ithvar(2) | bdd_ithvar(5));
    set1 = BDDTreeSet(Grid(3), 3, array<int>(3, 1,0,-1), enabled_cells);
    std::cout << "Set1 = ";
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;
    // test removal of variables
    indices = Vector<uint>(2);
    indices[0] = 0;
    indices[1] = 2;
    set2 = project_down(set1, indices);
    std::cout << "Set2 = ";
    for(BDDTreeSet::const_iterator it = set2.begin(); it != set2.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;
    
    BDDTreeSet set3(Grid(2), false);
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        Box bx = it->project(indices);
        set3.adjoin_lower_approximation(bx, 2);
    }
    std::cout << "Set3 = ";
    for(BDDTreeSet::const_iterator it = set3.begin(); it != set3.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;

    ARIADNE_TEST_EQUAL(set2, set3);
    
    // test reordering of variables
    indices = Vector<uint>(3);
    indices[0] = 2;
    indices[1] = 1;
    indices[2] = 0;
    set2 = project_down(set1, indices);
    std::cout << "Set2 = ";
    for(BDDTreeSet::const_iterator it = set2.begin(); it != set2.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;
    
    set3 = BDDTreeSet(Grid(3), false);
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        Box bx = it->project(indices);
        set3.adjoin_lower_approximation(bx, 2);
    }
    std::cout << "Set3 = ";
    for(BDDTreeSet::const_iterator it = set3.begin(); it != set3.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;

    ARIADNE_TEST_EQUAL(set2, set3);

    // test duplication of variables
    indices = Vector<uint>(3);
    indices[0] = 0;
    indices[1] = 1;
    indices[2] = 1;
    set2 = project_down(set1, indices);
    std::cout << "Set2 = ";
    for(BDDTreeSet::const_iterator it = set2.begin(); it != set2.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;
    
    set3 = BDDTreeSet(Grid(3), false);
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        Box bx = it->project(indices);
        set3.adjoin_lower_approximation(bx, 2);
    }
    std::cout << "Set3 = ";
    for(BDDTreeSet::const_iterator it = set3.begin(); it != set3.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;

    ARIADNE_TEST_EQUAL(set2, set3);

    // test projection of a minced set
    indices = Vector<uint>(3);
    indices[0] = 2;
    indices[1] = 0;
    indices[2] = 2;
    set1.mince(1);
    set2 = project_down(set1, indices);
    std::cout << "Set2 = ";
    for(BDDTreeSet::const_iterator it = set2.begin(); it != set2.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;
    
    set3 = BDDTreeSet(Grid(3), false);
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        Box bx = it->project(indices);
        set3.adjoin_lower_approximation(bx, 2);
    }
    set3.mince(1);
    std::cout << "Set3 = ";
    for(BDDTreeSet::const_iterator it = set3.begin(); it != set3.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;

    ARIADNE_TEST_EQUAL(set2, set3);
    
    // Test Luca's bug
    indices = Vector<uint>(1);
    indices[0] = 0;
    Grid grid(2);
    grid.set_length(0, 0.624);
    grid.set_length(1, 1.2);
    array<int> coordinates(2, 1,0);
    enabled_cells = bdd_nithvar(1) & ( (bdd_biimp(bdd_nithvar(0), bdd_ithvar(2)) & bdd_biimp(bdd_nithvar(3), bdd_ithvar(5)))
                        | (bdd_ithvar(0) & bdd_ithvar(2) & bdd_nithvar(3) & bdd_nithvar(4) & bdd_ithvar(5)) );
    set1 = BDDTreeSet(grid, 6, coordinates, enabled_cells);
    std::cout << "Set1 = " << set1 << std::endl;
    std::cout << "Set1 = ";
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        std::cout << *it << ", ";
    }
    
    set2 = project_down(set1, indices);
    std::cout << "Set2 = ";
    for(BDDTreeSet::const_iterator it = set2.begin(); it != set2.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;
    
    set3 = BDDTreeSet(Grid(1, 0.624), false);
    for(BDDTreeSet::const_iterator it = set1.begin(); it != set1.end(); it++) {
        Box bx = it->project(indices);
        set3.adjoin_lower_approximation(bx, 2);
    }
    std::cout << "Set3 = ";
    for(BDDTreeSet::const_iterator it = set3.begin(); it != set3.end(); it++) {
        std::cout << *it << ", ";
    }
    std::cout << std::endl;

    ARIADNE_TEST_EQUAL(set2, set3);

}

void test_restriction_difference() {
    ARIADNE_PRINT_TEST_COMMENT("Testing restriction with a ConstraintSet and a Checker.");
    //  Create the ImageSet to discretise
    RealVariable x("x");
    RealVariable y("y");
	List<RealVariable> varlist;
	varlist.append(x);
	varlist.append(y);
	RealExpression f_x = (y+1.0)*x + 5.0*y + 5.0;
	RealExpression f_y = -(y+1.0)*x + 5.0*y + 5.0;
	List<RealExpression> expr;
	expr.append(f_x);
	expr.append(f_y);
	VectorFunction f(expr,varlist);
	Box dom(2, -1.0,1.0, -1.0,1.0);
	ImageSet is1(dom, f);
    BDDTreeSet set1(2);
    set1.adjoin_outer_approximation(is1, 2);
    plot("test_bdd_set_restriction_is1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set1);
    
    // Create the ConstraintSet for the Restriction and the corresponding checker
    RealExpression invf_x = -x + 9.5;
	RealExpression invf_y = -y + 9.5;
    expr.clear();
    expr.append(invf_x);
	expr.append(invf_y);
	VectorFunction invf(expr,varlist);
	ConstraintSet cs1(invf, Box::upper_quadrant(2));
	ConstraintSetChecker checker1(cs1);

    // Testing outer restriction
    BDDTreeSet set0(0);
    ARIADNE_TEST_FAIL(set0.outer_restrict(cs1));
    ARIADNE_TEST_FAIL(set0.outer_restrict(checker1, 2));
    BDDTreeSet set3(3, true);
    ARIADNE_TEST_FAIL(set3.outer_restrict(cs1));    
    ARIADNE_TEST_FAIL(set0.outer_restrict(checker1, 2));
    BDDTreeSet set2(2, false); 
    set2.outer_restrict(cs1);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2.outer_restrict(checker1, 2);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2 = set1;
    set2.outer_restrict(cs1);
    plot("test_bdd_set_outer_restrict_cs1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set2);
    ARIADNE_TEST_ASSERT(set2.bounding_box().covers(Box(2, 0.0,9.5, 0.0,9.5)));
    ARIADNE_TEST_ASSERT(set2.subset(set1));
    ARIADNE_TEST_ASSERT(!set1.subset(set2));
    set3 = set1;
    set3.outer_restrict(checker1, 2);
    plot("test_bdd_set_outer_restrict_checker1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set3);
    ARIADNE_TEST_EQUAL(set2, set3);
    
    // Testing inner restriction
    ARIADNE_TEST_FAIL(set0.inner_restrict(cs1));
    ARIADNE_TEST_FAIL(set0.inner_restrict(checker1, 2));
    set3 = BDDTreeSet(3, true);
    ARIADNE_TEST_FAIL(set3.inner_restrict(cs1));    
    ARIADNE_TEST_FAIL(set3.inner_restrict(checker1, 2));    
    set2.clear(); 
    set2.inner_restrict(cs1);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2.inner_restrict(checker1, 2);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2 = set1;
    set2.inner_restrict(cs1);
    ARIADNE_TEST_ASSERT(!set2.empty());
    plot("test_bdd_set_inner_restrict_cs1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set2);
    ARIADNE_TEST_ASSERT(set2.bounding_box().inside(Box(2, -1.0,9.5, -1.0,9.5)));
    ARIADNE_TEST_ASSERT(set2.subset(set1));
    ARIADNE_TEST_ASSERT(!set1.subset(set2));
    set3 = set1;
    set3.inner_restrict(checker1, 2);
    plot("test_bdd_set_inner_restrict_checker1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set3);
    ARIADNE_TEST_EQUAL(set2, set3);

    ARIADNE_PRINT_TEST_COMMENT("Testing difference with a ConstraintSet and a Checker.");

    // Testing outer difference
    ARIADNE_TEST_FAIL(set0.outer_remove(cs1));
    ARIADNE_TEST_FAIL(set0.outer_remove(checker1, 2));
    set3 = BDDTreeSet(3, true);
    ARIADNE_TEST_FAIL(set3.outer_remove(cs1));    
    ARIADNE_TEST_FAIL(set0.outer_remove(checker1, 2));
    set2 = BDDTreeSet(2, false); 
    set2.outer_remove(cs1);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2.outer_remove(checker1, 2);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2 = set1;
    set2.outer_remove(cs1);
    plot("test_bdd_set_outer_remove_cs1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set2);
    BDDTreeSet set4(2);
    set4.adjoin_outer_approximation(set1.root_cell(), 2);
    set4.outer_restrict(cs1);
    set4 = difference(set1, set4);
    ARIADNE_TEST_EQUAL(set2,set4);
    set3 = set1;
    set3.outer_remove(checker1, 2);
    plot("test_bdd_set_outer_remove_checker1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set3);
    ARIADNE_TEST_EQUAL(set2, set3);
    
    // Testing inner restriction
    ARIADNE_TEST_FAIL(set0.inner_remove(cs1));
    ARIADNE_TEST_FAIL(set0.inner_remove(checker1, 2));
    set3 = BDDTreeSet(3, true);
    ARIADNE_TEST_FAIL(set3.inner_remove(cs1));    
    ARIADNE_TEST_FAIL(set3.inner_remove(checker1, 2));    
    set2.clear(); 
    set2.inner_remove(cs1);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2.inner_remove(checker1, 2);
    ARIADNE_TEST_ASSERT(set2.empty());
    set2 = set1;
    set2.inner_remove(cs1);
    ARIADNE_TEST_ASSERT(!set2.empty());
    plot("test_bdd_set_inner_remove_cs1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set2);
    set4.clear();
    set4.adjoin_outer_approximation(set1.root_cell(), 2);
    set4.inner_restrict(cs1);
    set4 = difference(set1, set4);
    ARIADNE_TEST_EQUAL(set2,set4);
    set3 = set1;
    set3.inner_remove(checker1, 2);
    plot("test_bdd_set_inner_remove_checker1",PlanarProjectionMap(2,0,1),set1.bounding_box(),Colour(1,0,1),set3);
    ARIADNE_TEST_EQUAL(set2, set3);

}


int main() {

    test_constructors();
    test_properties_subdivisions();
    test_predicates();
    test_operations();
    test_box_approximations();
    test_set_approximations();
    test_iterators_conversions_drawing();
    test_projection();
    test_restriction_difference();

    return ARIADNE_TEST_FAILURES;
}

