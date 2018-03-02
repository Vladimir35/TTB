#!/bin/bash
rake bin/ttbarxsec.exe
ttbarxsec.exe ../tt_PowhegP8.txt tt_PowhegP8.root --thread 1 --j 15 --J 600

