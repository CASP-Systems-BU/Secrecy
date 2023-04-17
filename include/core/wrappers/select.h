#ifndef SECRECY_SELECT_H
#define SECRECY_SELECT_H

#include "../../common/table.h"

//////////////////////////////////////////////////////////
///       Comparison        /      (Boolean - table)     /
////////////////////////////////////////////////////////// 
void select_EQ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                       const int &colTwo, const int &colOut);
void select_GT_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                       const int &colTwo, const int &colOut);
void select_GEQ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                   const int &colTwo, const int &colOut);
void select_EQZ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                   const int &colTwo, const int &colOut);
void select_GTZ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                   const int &colTwo, const int &colOut);
void select_GEQZ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                    const int &colTwo, const int &colOut);

//////////////////////////////////////////////////////////
///       Comparison        /      (Boolean - array)     /
//////////////////////////////////////////////////////////
ColumnData select_EQ_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_GT_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_GEQ_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_EQZ_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_GTZ_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_GEQZ_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);

//////////////////////////////////////////////////////////
///   Boolean Operations    /      (Boolean - table)     /
//////////////////////////////////////////////////////////
void
select_AND_table_b(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                   const int &colTwo, const int &colOut);

//////////////////////////////////////////////////////////
///   Boolean Operations    /      (Boolean - array)     /
//////////////////////////////////////////////////////////
ColumnData select_AND_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_OR_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_NOT_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows);

//////////////////////////////////////////////////////////
///  Arithmetic Operations  /    (Arithmetic - table)    /
////////////////////////////////////////////////////////// 
void select_ADD_table_a(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                   const int &colTwo, const int &colOut);
void select_SUB_table_a(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                   const int &colTwo, const int &colOut);
void select_MULT_table_a(const BShareTable &inputOne, const BShareTable &inputTwo, BShareTable &output, const int &colOne,
                    const int &colTwo, const int &colOut);

//////////////////////////////////////////////////////////
///  Arithmetic Operations  /    (Arithmetic - array)    /
//////////////////////////////////////////////////////////
ColumnData select_ADD_array_a(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_SUB_array_a(ColumnData inputOne, ColumnData inputTwo, const int &rows);
ColumnData select_MULT_array_a(ColumnData inputOne, ColumnData inputTwo, const int &rows);

#endif //SECRECY_SELECT_H
