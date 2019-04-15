/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 2009--2010  Nikolay Balov
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.gnu.org/licenses/gpl-2.0.html
 */

/*
 * thread.cpp
 *
 *  Created on: June 23, 2010
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

/* thread functionality */

#include "thread.h"
#include <stdio.h>

#if defined(WIN32) || defined(WIN64)

c_thread::c_thread() {
	m_hThread = NULL;
	InitializeCriticalSection(&m_thread_cs);
	m_thread_attr.nLength = sizeof(SECURITY_ATTRIBUTES);
	m_thread_attr.lpSecurityDescriptor = DEFAULT_SECURITY;
	m_thread_attr.bInheritHandle = FALSE;
	m_hStopThreadEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
	m_hJobDoneEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
}

c_thread::c_thread(THREAD_PROC ThreadProc, LPVOID lpParam) {

	m_hThread = NULL;
	InitializeCriticalSection(&m_thread_cs);
	m_thread_attr.nLength = sizeof(SECURITY_ATTRIBUTES);
	m_thread_attr.lpSecurityDescriptor = DEFAULT_SECURITY;
	m_thread_attr.bInheritHandle = FALSE;
	m_hStopThreadEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
	m_hJobDoneEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
	
	_start_thread(ThreadProc, lpParam);
}

c_thread::~c_thread() {
	_stop_thread();
	CloseHandle(m_hJobDoneEvent);
	CloseHandle(m_hStopThreadEvent);
	DeleteCriticalSection(&m_thread_cs);
}

void c_thread::_lock() {
	EnterCriticalSection(&m_thread_cs);
}

INT c_thread::_try_lock() {
	return TryEnterCriticalSection(&m_thread_cs);
}

void c_thread::_unlock() {
	LeaveCriticalSection(&m_thread_cs);
}

int c_thread::_wait_stop_event(int milliseconds) {
	int res;
	res = WaitForSingleObject(m_hStopThreadEvent, milliseconds);
	return res;
}

int c_thread::_join_thread() {
	int res;
	res = WaitForSingleObject(m_hJobDoneEvent, WAIT_INFINITE);
	return ERR_OK;
}

int c_thread::_start_thread(THREAD_PROC ThreadProc, LPVOID lpParam) {

	void *res = 0;
	UINT nId;
	if(m_hThread && WaitForSingleObject(m_hThread, 0) != WAIT_OBJECT_0) {
		return ERR_THREAD;
	}
	ResetEvent(m_hStopThreadEvent);
	ResetEvent(m_hJobDoneEvent);
	__TRY_FINALLY {
		_lock();
		m_hThread = (HANDLE)_beginthreadex(NULL, 0, ThreadProc, lpParam, 0, &nId);
		res = m_hThread;
	}
	__FINALLY {
		_unlock();
	}
	__END_FINALLY
	
	return ((res != 0) ? ERR_OK : ERR_THREAD);
}

int c_thread::_stop_thread() {

	if(!m_hThread)
		return ERR_OK;

	SetEvent(m_hStopThreadEvent);
	WaitForSingleObject(m_hThread, INFINITE);
	
	__TRY_FINALLY {
		_lock();
		CloseHandle((HANDLE)m_hThread);
		m_hThread = NULL;
	}
	__FINALLY {
		_unlock();
	}
	__END_FINALLY

	return ERR_OK;
}

int c_thread::_exit_thread(void* pexitcode) {
	SetEvent(m_hJobDoneEvent);
	return ERR_OK;
}

/*  
 THREAD_PROC_DEFINE(SampleThreadProc, pParam){

 c_thread * pThreadObj = (c_thread*)pParam;
 while(WaitForSingleObject(pThreadObj->m_hStopThreadEvent, 0) == WAIT_TIMEOUT){
 __TRY_FINALLY{
 // do some work ...
 }
 __FINALLY{
 pThreadObj->_unlock();
 }
 __END_FINALLY
 }        
 return exit_code;
 }
 */

#else //#ifdef _POSIX_THREADS

c_thread::c_thread() {
	m_hThread = 0;
	pthread_mutex_init(&m_thread_cs, NULL);
	pthread_attr_init(&m_thread_attr);
	pthread_attr_setdetachstate(&m_thread_attr, PTHREAD_CREATE_JOINABLE);

	pthread_cond_init(&m_stop_event, NULL);
	pthread_mutex_init(&m_stop_mutex, NULL);
}

c_thread::c_thread(THREAD_PROC ThreadProc, LPVOID lpParam) {
	m_hThread = 0;
	pthread_mutex_init(&m_thread_cs, NULL);
	pthread_attr_init(&m_thread_attr);
	pthread_attr_setdetachstate(&m_thread_attr, PTHREAD_CREATE_JOINABLE);

	pthread_cond_init(&m_stop_event, NULL);
	pthread_mutex_init(&m_stop_mutex, NULL);

	_start_thread(ThreadProc, lpParam);
}

c_thread::~c_thread() {
	_stop_thread();
	pthread_attr_destroy(&m_thread_attr);
	pthread_cond_destroy(&m_stop_event);
	pthread_mutex_destroy(&m_stop_mutex);
	pthread_mutex_destroy(&m_thread_cs);
}

void c_thread::_lock() {
	pthread_mutex_lock(&m_thread_cs);
}

int c_thread::_try_lock() {
	return pthread_mutex_trylock(&m_thread_cs);
}

void c_thread::_unlock() {
	pthread_mutex_unlock(&m_thread_cs);
}

int c_thread::_wait_stop_event(int milliseconds) {
	int res;
	if (milliseconds == WAIT_INFINITE)
		milliseconds = 0x7fffffff;
	timespec abstime;
	//clock_gettime(CLOCK_REALTIME, &abstime);
	abstime.tv_sec = (int) (milliseconds / 1000);
		milliseconds -= 1000 * (int) (milliseconds / 1000);
	abstime.tv_nsec = milliseconds * 1000000;
	pthread_mutex_lock(&m_stop_mutex);
	res = pthread_cond_timedwait(&m_stop_event, &m_stop_mutex, &abstime);
	pthread_mutex_unlock(&m_stop_mutex);
	return res;
}

int c_thread::_join_thread() {
	void *pstatus = 0;
	int res = 0;
	_lock();
	if (!m_hThread) {
		_unlock();
		return res;
	}
	_unlock();
	res = pthread_join(m_hThread, &pstatus);
	if(!res) {
		_lock();
		m_hThread = 0;
		_unlock();
	}
	return res;
}

int c_thread::_start_thread(THREAD_PROC ThreadProc, LPVOID lpParam) {
	int res;
	_lock();
	if (m_hThread) {
		_unlock();
		return ERR_OK;
	}
	res = pthread_create(&m_hThread, &m_thread_attr, ThreadProc, lpParam);
	_unlock();
	return ((res != 0) ? ERR_OK : ERR_THREAD);
}

int c_thread::_stop_thread() {
	int res = ERR_OK;
	_lock();
	if (!m_hThread) {
		_unlock();
		return res;
	}
	_unlock();
	pthread_cond_signal(&m_stop_event);
	res = _join_thread();
	if(res) {
		res = pthread_cancel(m_hThread);
		_lock();
		m_hThread = 0;
		_unlock();
	}
	if(!res)
		return ERR_OK;
	return ERR_STOP_THREAD;
}

int c_thread::_exit_thread(void* pexitcode) {
	pthread_exit(pexitcode);
	return ERR_OK;
}

/*  
 THREAD_PROC_DEFINE(SampleThreadProc, pParam){

 c_thread * pThreadObj = (c_thread*)pParam;
 while(1){
 pthread_mutex_lock(&pThreadObj->m_stop_mutex);
 if(pThreadObj->_wait_event(pThreadObj->m_stop_event, SOME_TIME) == 0)
 break;   
 // do some work ...
 pthread_mutex_unlock(&pThreadObj->m_stop_mutex);
 }        
 THREAD_EXIT(exit_code);
 }
 */

#endif /* _POSIX_THREADS */

/* end of file */
