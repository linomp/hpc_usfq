/*
 * auxFunctions.h
 *
 *  Created on: May 4, 2016
 *      Author: reid
 */

#ifndef AUXFUNCTIONS_H_
#define AUXFUNCTIONS_H_

#include "definitions.h"


Vector4d qfromR(const Matrix3d & R) {
	Vector4d _quat;
	double tr = R.trace();
	_quat(3) = sqrt(tr+1.)/2.;
	_quat(0) = (R(0,2)-R(2,0))/(4*_quat(3));
	_quat(1) = (R(1,2)-R(2,1))/(4*_quat(3));
	_quat(2) = (R(0,1)-R(1,0))/(4*_quat(3));
	return _quat;
}

Matrix3d Rfromq(const Vector4d & quat) {
	Matrix3d _rotMatrix;
	_rotMatrix(0,0) = -quat(0)*quat(0) + quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
	_rotMatrix(0,1) = -2*(quat(0)*quat(1) - quat(2)*quat(3));
	_rotMatrix(0,2) =  2*(quat(1)*quat(2) + quat(0)*quat(3));
	_rotMatrix(1,0) = -2*(quat(0)*quat(1) + quat(2)*quat(3));
	_rotMatrix(1,1) =  quat(0)*quat(0) - quat(1)*quat(1) - quat(2)*quat(2) + quat(3)*quat(3);
	_rotMatrix(1,2) = -2*(quat(0)*quat(2) - quat(1)*quat(3));
	_rotMatrix(2,0) =  2*(quat(1)*quat(2) - quat(0)*quat(3));
	_rotMatrix(2,1) = -2*(quat(0)*quat(2) + quat(1)*quat(3));
	_rotMatrix(2,2) = -quat(0)*quat(0) - quat(1)*quat(1) + quat(2)*quat(2) + quat(3)*quat(3);
	return _rotMatrix;
}



#endif /* AUXFUNCTIONS_H_ */
