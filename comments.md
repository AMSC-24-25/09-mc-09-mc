# comments

You have not provided a real report, just an extensive `README.md` file. This can be enough, but in fact there is something missing: an explanation on how to compile the code and run some tests. Luckily, I understood what to do.

I would have been nicer if some of the examples could read the parameters from a file, so that one cha chenge them without recompiling the code.

As for the rest. The code looks well written. I have put some notes marked `@note` in the sources. 

Like most of students of this year, you do not like `hpp` for C++ header files. `h` is indeed fine but it is nice for the reader to understand that it is a C++ header file. 

Remember that a more professional code should be organised as a library, and when you rune make install you should install the library in the specified system directories. Typically, in a development code, you put the header files in /usr/local/include and the precompiled libraries in /usr/local/lib (if you have them). But this is not mandatory, and you can also install them in another place.

Sometimes the const qualifier for methods is missing. In your case is not a big issue, however remeber that a user expect to be able to call on a contant object any method that does not change the state of the object. 

Fore the rest, the code is well written with a good use of some modern c++ features. I have seen a lot of commits, so I guess you have worked a lot on this project.


