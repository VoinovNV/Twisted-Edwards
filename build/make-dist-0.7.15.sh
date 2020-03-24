#/bin/bash
mkdir -p libakrypt-0.7.15/source
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/source/* libakrypt-0.7.15/source
mkdir -p libakrypt-0.7.15/examples
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/examples/example-*.c libakrypt-0.7.15/examples
mkdir -p libakrypt-0.7.15/tests
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/tests/test-*.c libakrypt-0.7.15/tests
mkdir -p libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/00-introduction.md libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/01-install-guide.md libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/02-principles.md libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/asn1.md libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/Doxyfile.in libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/refman_footer.in libakrypt-0.7.15/doc
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/doc/refman_header.in libakrypt-0.7.15/doc
mkdir -p libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/DetectFlags.cmake libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/DetectHeaders.cmake libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/DetectIntrinsics.cmake libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/DetectLibraries.cmake libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/MakeDist.cmake libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/MakeDoc.cmake libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/cmake/MakeInstall.cmake libakrypt-0.7.15/cmake
mkdir -p libakrypt-0.7.15/akrypt
mkdir -p libakrypt-0.7.15/cmake
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/CMakeLists.txt libakrypt-0.7.15
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/Readme.md libakrypt-0.7.15
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/Changelog.md libakrypt-0.7.15
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/libakrypt.conf libakrypt-0.7.15
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/libakrypt.rc.in libakrypt-0.7.15
cp -fL --preserve=all /home/user/Projects/Twisted_Edwards/libakrypt-0.x/LICENSE libakrypt-0.7.15
tar -cjvf libakrypt-0.7.15.tar.bz2 libakrypt-0.7.15
rm -R libakrypt-0.7.15
akrypt i -r -o libakrypt-0.7.15.streebog256 libakrypt-0.7.15.tar.bz2
