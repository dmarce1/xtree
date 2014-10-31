/*
 * new.cpp
 *
 *  Created on: Oct 28, 2014
 *      Author: dmarce1
 */

#include <new>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>

#include "tbb/cache_aligned_allocator.h"
#include "tbb/scalable_allocator.h"

static tbb::scalable_allocator<char> small_alloc;
static tbb::cache_aligned_allocator<char> big_alloc;
static const std::size_t cache_line_size = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);

static void* allocate(std::size_t n) {
	if (n >= cache_line_size) {
		return big_alloc.allocate(n);
	} else {
		return small_alloc.allocate(n);
	}
}

static int deallocate(void* ptr) {
	free(ptr);
}

void* operator new(std::size_t n) {
	return allocate(n);
}

void* operator new[](std::size_t n) {
	return allocate(n);

}

void* operator new(std::size_t n, const std::nothrow_t& nothrow_value) noexcept
{
	void *ptr;
	try {
		ptr = allocate(n);
	} catch (...) {
		ptr = nullptr;
	}
}

void* operator new[](std::size_t n, const std::nothrow_t& nothrow_value) noexcept
{
	void *ptr;
	try {
		ptr = allocate(n);
	} catch (...) {
		ptr = nullptr;
	}

}

void operator delete(void* ptr) {
	deallocate(ptr);

}

void operator delete[](void* ptr) {
	deallocate(ptr);
}

void operator delete(void* ptr, const std::nothrow_t& ) noexcept {
	deallocate(ptr);

}

void operator delete[](void* ptr, const std::nothrow_t& ) noexcept {
	deallocate(ptr);
}
