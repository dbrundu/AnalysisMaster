/*----------------------------------------------------------------------------
 *
 *   Copyright (C) 2018 D.Brundu, A.Contu et al.
 *
 *   This file is part of D0->hhmumu analysis package.
 *   It uses Hydra as external header-only software.
 * 
 *---------------------------------------------------------------------------*/

/*
 * utils.h
 *
 *  Created on: 30/04/2018
 *      Author: augalves
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>

#define VERBOSE_LINE(flag, message)\
if(flag){\
	std::cout<< "\033[1;35mVerbose: \033[0m"\
			 << message << std::endl;\
}


#define INFO_LINE(message)\
	std::cout<< "\033[1;34mInfo: \033[0m"\
			 << message << std::endl;\

#define WARNING_LINE(message)\
	std::cout<< "\033[1;34mWarning: \033[0m"\
			 << message << std::endl;\




#endif /* UTILS_H_ */
