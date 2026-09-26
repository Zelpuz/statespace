Examples
========

The illustrations of DK ch. 8 and the two models of DK §3.10, as Fortran
programs in ``example/``. Each page compares the results with the values DK
print. Run them from the repository root after fetching the data::

   python data/fetch_dk_data.py
   fpm run --profile release --example dk_8_2_seatbelt

Where DK's printed values could be traced, the programs reproduce them. The
pages also record misprints and conventions that differ from ours.

.. toctree::
   :maxdepth: 1

   dk_8_2
   dk_8_3
   dk_8_4
   dk_8_5
   dk_8_6
   dk_3_10_2
   dk_3_10_3
