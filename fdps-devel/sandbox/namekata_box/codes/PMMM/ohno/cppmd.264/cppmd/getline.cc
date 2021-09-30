// Copyright (c) 2010 RIKEN. All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//    1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//    2. Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//    3. The name of the author may not be used to endorse or promote
//    products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY RIKEN ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif  // HAVE_CONFIG_H
#include "getline.h"

#ifndef HAVE_GETLINE

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {

const int kINITIAL_BUFFER_SIZE = 127;

}  // namespace

ssize_t getline(char **lineptr, size_t *n, FILE *stream) {
  return getdelim(lineptr, n, '\n', stream);
}

ssize_t getdelim(char **lineptr, size_t *n, int delim, FILE *stream) {
  // validate arguments
  if (n == NULL || lineptr == NULL || stream == NULL) {
    errno = EINVAL;
    return -1;
  }

  // update *lineptr and *n as necessary
  if (*lineptr == NULL
      || *n < static_cast<size_t>(kINITIAL_BUFFER_SIZE)) {
    size_t newsize = kINITIAL_BUFFER_SIZE;
    char *newptr = static_cast<char*>(realloc(*lineptr, newsize));
    if (newptr == NULL) return -1;
    *lineptr = newptr;
    *n = newsize;
  }

  ssize_t readlen = 0;
  (*lineptr)[readlen] = '\0';
  for (;;) {
    if (feof(stream)) return -1;
    // update *lineptr and *n as necessary
    if (*n - readlen <= 1) {
      size_t newsize = *n * 2 + 1;
      if (newsize < *n) return -1;
      char *newptr = static_cast<char*>(realloc(*lineptr, newsize));
      if (newptr == NULL) return -1;
      *lineptr = newptr;
      *n = newsize;
    }
    char *s = *lineptr + readlen;
    while (*n - readlen - 1 > 0) {
      int ch = getc(stream);
      if (ch != EOF) {
        *s++ = static_cast<char>(ch);
        ++readlen;
      }
      if (ch == delim || ch == EOF) {
        if (ferror(stream)) return -1;
        *s = '\0';
        break;
      }
    }
    if ((*lineptr)[readlen - 1] == delim || feof(stream))
      break;
  }
  return readlen;
}

#endif  // HAVE_GETLINE
