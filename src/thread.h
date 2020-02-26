/*
 *  catnet : categorical Bayesian network inference
 *  Copyright (C) 1999--2010  Nikolay Balov
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
 * thread.h
 *
 *  Created on: June 23, 2010
 *      Author: nbalov
 */

/* 
 * version 1.15.1  12dec2016
 */

/* thread functionality */

#ifndef _THREAD_H_
#define _THREAD_H_

#if defined(WIN32) || defined(WIN64)

#define _OBJC_NO_COM 
#define NOGDI
#include <windows.h>
#include <process.h>

typedef HANDLE THREAD;
typedef CRITICAL_SECTION MUTEX;
typedef HANDLE EVENT;

typedef unsigned (__stdcall * THREAD_PROC)(void*);
#define THREAD_PROC_DEFINE(proc_name, proc_param) unsigned __stdcall proc_name(void *proc_param) 

#define CREATE_EVENT(event)			(event = CreateEvent(NULL, TRUE, FALSE, NULL))
#define DESTROY_EVENT(event)			CloseHandle(event)
#define SET_EVENT(event)			SetEvent(event)
#define RESET_EVENT(event)			ResetEvent(event)
#define WAIT_EVENT(event, mutex, time)		WaitForSingleObject(&event, &time)

#define MUTEX_LOCK(pmutex)			EnterCriticalSection((CRITICAL_SECTION*)pmutex)
#define MUTEX_UNLOCK(pmutex)			LeaveCriticalSection((CRITICAL_SECTION*)pmutex)

#define MUTEX_INIT(mutex)			InitializeCriticalSection((CRITICAL_SECTION*)&mutex)
#define MUTEX_DESTROY(mutex)			DeleteCriticalSection((CRITICAL_SECTION*)&mutex)

#define THREAD_EXIT(status)			return status;

#define WAIT_INFINITE       INFINITE

#define DEFAULT_SECURITY	NULL

#define __TRY_EXCEPT		{
#define __EXCEPT(expr)		}if(expr)	\
							{
#define __END_EXCEPT		}

#define __TRY_FINALLY		{
#define __FINALLY			}	\
							{
#define __END_FINALLY		}
#define __LEAVE_FINALLY


#else //defined(_POSIX_THREADS)

#include <pthread.h>

typedef pthread_t		THREAD;
typedef pthread_mutex_t		MUTEX;
typedef pthread_cond_t		EVENT;

typedef void * (*THREAD_PROC)(void*);
#define THREAD_PROC_DEFINE(proc_name, proc_param) void* proc_name(void* proc_param)

/* static initializing cond = PTHREAD_COND_INITIALIZER */
#define CREATE_EVENT(event)			pthread_cond_init(&event, NULL)
#define DESTROY_EVENT(event)			pthread_cond_destroy(&event)
#define SET_EVENT(event)			pthread_cond_signal(&event)
#define RESET_EVENT(event)        
#define BROADCAST_EVENT(even)			pthread_cond_broadcast(&event)
#define WAIT_EVENT(event, mutex, time)		pthread_cond_timedwait(&event, &mutex, &time)

#define MUTEX_LOCK(pmutex)			pthread_mutex_lock(pmutex)
#define MUTEX_UNLOCK(pmutex)			pthread_mutex_unlock(pmutex)

#define MUTEX_INIT(mutex)			pthread_mutex_init(&mutex, NULL)
#define MUTEX_DESTROY(mutex)			pthread_mutex_destroy(&mutex)

#define WAIT_INFINITE   -1
#define WAIT_TIMEOUT    ETIMEDOUT

#define THREAD_EXIT(status)              pthread_exit(status)

#define __TRY_EXCEPT		{
#define __EXCEPT(expr)		}if(expr)	\
							{
#define __END_EXCEPT		}

#define __TRY_FINALLY		{
#define __FINALLY			}	\
							{
#define __END_FINALLY		}
#define __LEAVE_FINALLY

#endif /* _POSIX_THREADS */

#define ERR_OK			0
#define ERR_THREAD		-1
#define ERR_START_THREAD	-2
#define ERR_STOP_THREAD		-3
#define LPVOID			void*

class c_thread {

	/* attributes */
protected:

#if defined(WIN32) || defined(WIN64)
	
	SECURITY_ATTRIBUTES	m_thread_attr;
	HANDLE			m_hThread;
	CRITICAL_SECTION	m_thread_cs;
	HANDLE			m_hStopThreadEvent;
	HANDLE			m_hJobDoneEvent;
#else //ifdef _POSIX_THREADS
	pthread_attr_t		m_thread_attr;
	pthread_t		m_hThread;
	pthread_mutexattr_t	m_mutex_attr;
	pthread_mutex_t		m_thread_cs;
	pthread_cond_t		m_stop_event;
	pthread_mutex_t		m_stop_mutex;
#endif

	/* members */
public:
	c_thread();
	c_thread(THREAD_PROC ThreadProc, LPVOID lpParam);
	virtual ~c_thread();

	/* m_thread_sc control */
	void _lock();
	int _try_lock();
	void _unlock();

	int _wait_stop_event(int milliseconds);
	int _join_thread();
	int _exit_thread(void* pexitcode);
	int _start_thread(THREAD_PROC ThreadProc, LPVOID lpParam);
	int _stop_thread();

	int _is_running() {
		return (m_hThread != 0);
	}
};

#endif	/* _THREAD_H_ */
