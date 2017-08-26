Installation
============

.. note:: Installation via `pip` will be available only once the paper is published.

Using :py:obj:`pip`:

.. code-block:: bash

   pip install planetplanet

To upgrade:

.. code-block:: bash

   pip install -U --no-deps planetplanet

From source:

.. code-block:: bash

  git clone git@github.com:rodluger/planetplanet.git
  cd planetplanet
  git submodule init && git submodule update
  python setup.py develop