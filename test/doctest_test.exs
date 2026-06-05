defmodule Orbis.DoctestTest do
  use ExUnit.Case

  doctest Orbis
  doctest Orbis.RF
  doctest Orbis.GNSS.Constellation
  doctest Orbis.Format.TLE
  doctest Orbis.Format.OMM
  doctest Orbis.GNSS.Signal.CA
  doctest Orbis.GNSS.Navigation.LNAV
end
