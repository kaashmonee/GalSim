Changes from v2.2 to v2.3
=========================


Deprecated Features
-------------------


API Changes
-----------


Config Updates
--------------

- Fixed a few issues with parsing truth catalog items and, in general, places
  where a Current item might be used before it was parsed. (#1083)
- Added value-type-specific type names like Random_float, Random_Angle, etc.
  to help in some cases where the config processing cannot know what value
  type to use for further processing.  (#1083)
- Fixed a subtle issue in Eval string processing if one Current item is a
  shorter version of another one.  e.g. @gal and @gal.index.  It had been
  the case that the longer one might not be parsed properly. (#1083)


New Features
------------


Bug Fixes
---------

- Fixed horner and horner2d when inputs are complex. (#1054)
