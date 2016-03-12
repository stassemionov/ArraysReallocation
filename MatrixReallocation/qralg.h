#ifndef _QRALG_H_
#define _QRALG_H_

#include "taskdata.h"

// ����������� QR-�������� � WY-�����������.
double* QR_WY_standard(double* A, const int N, const int r);

// QR-�������� � WY-�����������.
// ��������� ������ ������������ � ���������.
double* QR_WY_tiled(double* A, const TaskData& parameters);

// QR-�������� � WY-�����������.
// ������� � �������������� ����������� ������.
double* QR_WY_block(double* A, const TaskData& parameters);

#endif  // _QRALG_H_