 
OCVWarp v2.50
-----------------
Bugfixes. Output more closely matches GL_warp for transforms 4 and 5. The Windows builds are via Appveyor.


 OCVWarp v2.11
-----------------
Replacing the auto-generated AppImage which was inadvertently deleted in the previous release. Travis still doesn't include the tag of the release in the AppImage name, only the commit hash.


 OCVWarp v2.10
-----------------
This release fixes a couple of bugs - the seam at longitude=0 and angley not working in transformtype=0.

 
 OCVWarp v2.00
-----------------
This release adds support for transformtype=5 - Equirectangular 360 video to warped
Also added feature to ask for output filename, allowing output file path and container (format) type to be changed.

The wiki has been updated with basic info related to usage etc.

 OCVWarp v1.45 - Windows release
 -------------------------------

Windows exe with special characters suppressed on terminal. Tested working on Windows 10 64-bit and Windows 7 32-bit. Compiled with mingw32.

Change the ini file parameters to suit your needs.

Transformtype.txt gives details of transform types supported.

 OCVWarp v1.42 - Windows release
--------------------------------

This release includes a working Windows binary. Tested OK on Windows 10, not running on Windows 7. This has support for transformtype=4.

More details in transformtype.txt.

Download the zip file into a single directory and run the OCVWarp-1.41-x86Win.exe file.


 OCVWarp v1.40
------------------

This release includes support for transformtype=4. A Linux AppImage - pre-built static binary for recent x64 Linux distros - is included. See transformtype.txt for details.

The getfourcc binary is a helper, lists possible output codec fourcc values. Its output is listed in the appimagefourcc.txt file. See previous release for getfourcc binary.

Download the ini file and the OCVWarp Applmage, make the AppImage executable, and run it.

chmod +x OCVWarp-1.40-x86_64.AppImage


 OCVWarp v1.30
---------------

This release includes AppImages - pre-built static binaries for recent x64 Linux distros. Includes working functionality for converting Equirectangular to Fisheye and vice versa. Also added ability to choose output codec. Ini file is changed accordingly.

The getfourcc binary is a helper, lists possible output codec fourcc values. Its output is listed in the appimagefourcc.txt file.

Download the ini file and the OCVWarp Applmage, make the AppImage executable, and run it.

chmod +x OCVWarp1.30-x86_64.AppImage

 OCVWarp v1.20
----------------

This release includes an AppImage - pre-built static binary for recent x64 Linux distros. Includes working functionality for converting Equirectangular to Fisheye and vice versa.

 OCVWarpNorth v1.02
----------------------

This release has Windows binaries for OCVWarpNorth included. Please unzip the file, install the VC_redist_x64 files if needed. Built on appveyor.com.
