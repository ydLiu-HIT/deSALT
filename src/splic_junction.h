/*************************************************************************
	> File Name: splic_junction.h
	> Author: 
	> Mail: 
	> Created Time: 2017年11月17日 星期五 22时15分06秒
 ************************************************************************/

#ifndef _SPLIC_JUNCTION_H
#define _SPLIC_JUNCTION_H

void acceptor_signals_detected(uint8_t* ref, uint16_t length, uint32_t detected_start, int16_t *s1, int16_t *s2, uint8_t strand);

void donor_signals_detected(uint8_t* ref, uint16_t length, uint32_t detected_start, int16_t *s1, int16_t *s2, uint8_t strand);
#endif
