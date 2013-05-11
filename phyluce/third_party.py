import os

def which(name, flags=os.X_OK):
    """
    Copyright (c) 2001-2008
    Allen Short
    Andrew Bennetts
    Apple Computer, Inc.
    Benjamin Bruheim
    Bob Ippolito
    Canonical Limited
    Christopher Armstrong
    David Reid
    Donovan Preston
    Eric Mangold
    Itamar Shtull-Trauring
    James Knight
    Jason A. Mobarak
    Jean-Paul Calderone
    Jonathan Lange
    Jonathan D. Simms
    Jargen Hermann
    Kevin Turner
    Mary Gardiner
    Matthew Lefkowitz
    Massachusetts Institute of Technology
    Moshe Zadka
    Paul Swartz
    Pavel Pergamenshchik
    Ralph Meijer
    Sean Riley
    Software Freedom Conservancy
    Travis B. Hartwell
    Thomas Herve
    Eyal Lotem
    Antoine Pitrou
    Andy Gayton

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    Search PATH for executable files with the given name.
    
    On newer versions of MS-Windows, the PATHEXT environment variable will be
    set to the list of file extensions for files considered executable. This
    will normally include things like ".EXE". This fuction will also find files
    with the given name ending with any of these extensions.

    On MS-Windows the only flag that has any meaning is os.F_OK. Any other
    flags will be ignored.
    
    @type name: C{str}
    @param name: The name for which to search.
    
    @type flags: C{int}
    @param flags: Arguments to L{os.access}.
    
    @rtype: C{list}
    @param: A list of the full paths to files found, in the
    order in which they were found.
    """
    result = []
    exts = filter(None, os.environ.get('PATHEXT', '').split(os.pathsep))
    path = os.environ.get('PATH', None)
    if path is None:
        #return []
        raise EnvironmentError("Cannot find $PATH variable.")
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if os.access(p, flags):
            result.append(p)
        for e in exts:
            pext = p + e
            if os.access(pext, flags):
                result.append(pext)
    if result == []:
        raise EnvironmentError("Cannot find {0}.  Make sure {0} is in your $PATH.".format(name))
    return result
