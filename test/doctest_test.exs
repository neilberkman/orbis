defmodule Orbis.DoctestTest do
  use ExUnit.Case

  doctest Orbis
  doctest Orbis.RF
  doctest Orbis.Format.TLE
  doctest Orbis.Format.OMM
end
