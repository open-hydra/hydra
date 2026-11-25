# ---------------------------------------------------------------------------
#  CTest definitions of hydra: one ctest per case of test/<application>/, for the application this build makes
#  (EXE_NAME). Every test is registered with hydra_test() below: a name, a one-line DESC, the KIND of check and its
#  LABELS. Its tier is not written here: it comes from test/tiers.txt (one line per test), because the tier follows
#  from the measured time of the test, not from what it checks:
#
#    light    under 5 s     pre-push (ctest -L light)
#    medium   5 s to 100 s  every PR (ctest -L "light|medium")
#    heavy    100 s and up  local, on request (ctest -L heavy)
#
#  Labels of a test: its tier; its kind (exact identity contract robustness reference); its application (MI2).
#  The runs go to <build>/test/<application>/; each test's PASS/FAIL lines to <that folder>/<case>.check.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
#  Tiers: test/tiers.txt, "<Name> <light|medium|heavy> <seconds>" per line (the seconds at one thread). A change of
#  the file reconfigures.
# ---------------------------------------------------------------------------
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tiers.txt")
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/tiers.txt" _tier_lines)
set(_tier_names "")
foreach(_line IN LISTS _tier_lines)
    if(_line MATCHES "^[ \t]*(#|$)")
        continue()
    endif()
    if(NOT _line MATCHES "^([^ \t]+)[ \t]+(light|medium|heavy)[ \t]+([0-9.]+|-)[ \t]*$")
        message(FATAL_ERROR "test/tiers.txt: cannot read the line '${_line}'")
    endif()
    set(HYDRA_TIER_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    list(APPEND _tier_names ${CMAKE_MATCH_1})
endforeach()

# ---------------------------------------------------------------------------
#  hydra_test(<Name>
#      DESC "<what it checks, one line, no semicolon>"
#      KIND <exact|identity|contract|robustness|reference>
#      LABELS <labels>
#      COMMAND <program> <arguments>
#      [TIMEOUT <s>] [SETUP <fixture>] [REQUIRES <fixture>...] [SERIAL])
#  SETUP: the test's run is the base another test reads (ctest runs it first); REQUIRES: the bases it reads;
#  SERIAL: no other test runs beside it.
# ---------------------------------------------------------------------------
function(hydra_test NAME)
    cmake_parse_arguments(T "SERIAL" "KIND;DESC;TIMEOUT;SETUP" "LABELS;COMMAND;REQUIRES" ${ARGN})
    if(NOT T_DESC OR NOT T_KIND OR NOT T_LABELS OR NOT T_COMMAND)
        message(FATAL_ERROR "hydra_test(${NAME}): DESC, KIND, LABELS and COMMAND are required")
    endif()
    if(NOT T_KIND MATCHES "^(exact|identity|contract|robustness|reference)$")
        message(FATAL_ERROR "hydra_test(${NAME}): unknown KIND '${T_KIND}'")
    endif()
    if(NOT DEFINED HYDRA_TIER_${NAME})
        message(FATAL_ERROR "hydra_test(${NAME}): no line for it in test/tiers.txt")
    endif()
    set(_tier ${HYDRA_TIER_${NAME}})
    if(NOT T_TIMEOUT)
        if(_tier STREQUAL "light")
            set(T_TIMEOUT 60)
        elseif(_tier STREQUAL "medium")
            set(T_TIMEOUT 600)
        else()
            set(T_TIMEOUT 7200)
        endif()
    endif()
    add_test(NAME ${NAME} COMMAND ${T_COMMAND})
    set_tests_properties(${NAME} PROPERTIES LABELS "${_tier};${T_KIND};${T_LABELS}" TIMEOUT ${T_TIMEOUT}
                         ENVIRONMENT "KMP_STACKSIZE=100M")
    if(T_SETUP)
        set_tests_properties(${NAME} PROPERTIES FIXTURES_SETUP ${T_SETUP})
    endif()
    if(T_REQUIRES)
        set_tests_properties(${NAME} PROPERTIES FIXTURES_REQUIRED "${T_REQUIRES}")
    endif()
    if(T_SERIAL)
        set_tests_properties(${NAME} PROPERTIES RUN_SERIAL TRUE)
    endif()
    set_property(GLOBAL APPEND PROPERTY HYDRA_TEST_NAMES ${NAME})
endfunction()

# ---------------------------------------------------------------------------
#  hydra-MI2: the cases of test/MI2/cases (test/MI2/run.sh, one case per test) and the validation of
#  test/MI2/validation (each its own script)
# ---------------------------------------------------------------------------
if(EXE_NAME STREQUAL "MI2")
    set(_exe "${CMAKE_SOURCE_DIR}/bin/hydra_MI2")
    set(_runs "${CMAKE_CURRENT_BINARY_DIR}/MI2")
    set(_val "${CMAKE_CURRENT_SOURCE_DIR}/MI2/validation")
    # mi2_case(<case> <KIND> "<DESC>" [hydra_test options]): MI2.<case>, run.sh on that case
    function(mi2_case CASE KIND DESC)
        hydra_test(MI2.${CASE} KIND ${KIND} DESC "${DESC}" LABELS MI2 ${ARGN}
            COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/MI2/run.sh" -x "${_exe}" -s "${_runs}" -t 1 ${CASE})
    endfunction()

    mi2_case(igloo-evap          exact    "IGLOO drops evaporate into the gas: mass, species and energy books")
    mi2_case(igloo-evap-sa       exact    "igloo-evap with a RANS gas: the books and the turbulence variable left alone")
    mi2_case(igloo-two-materials exact    "two IGLOO materials each into its own species")
    mi2_case(igloo-burn          exact    "IGLOO drops burn into the gas against closed forms")
    mi2_case(igloo-al            exact    "IGLOO aluminium carries its oxide and generates alumina in ICE"
             SETUP MI2.igloo-al)
    mi2_case(igloo-al-release    exact    "igloo-al with the carried oxide released into ICE at burnout")
    mi2_case(ice-exchange        exact    "ICE drag and heat to the gas against closed forms")
    mi2_case(ice-burn            exact    "ICE drops burn into the gas against closed forms")
    mi2_case(ice-burnout         exact    "ICE drops burn out: the books of the removal")
    mi2_case(ice-burnout-steady  identity "ice-burnout in pseudo-time: the converged field independent of ICE's step")
    mi2_case(ice-al              exact    "ICE aluminium generates alumina in another family: one step against closed forms"
             SETUP MI2.ice-al)
    mi2_case(ice-igloo-burn      exact    "ICE and IGLOO burning in one gas: the books per phase")
    mi2_case(ch-gas-ice          exact    "channel gas to ICE: ICE's gain and the gas's rate against ice-al"
             REQUIRES MI2.ice-al)
    mi2_case(ch-gas-igloo        exact    "channel gas to IGLOO: the lag and the books against igloo-al"
             REQUIRES MI2.igloo-al)
    mi2_case(ch-steady           exact    "channels IGLOO-ICE and gas-ICE in pseudo-time: the rates")
    mi2_case(refusals            contract "[MI2-Channels] lines refused at setup")
    mi2_case(twins               identity "one burning section in ICE and IGLOO: energy per kg and law constant agree")
    mi2_case(burnout             exact    "ice-burnout variants: long step, frozen, cold; d-stop beside d-burnout refused")
    mi2_case(routing             contract "phase routing: which solver takes which phase file")
    mi2_case(ini-reread          identity "input.ini re-read during the run changes nothing")
    mi2_case(floor-guard         contract "the consumed-species floor guard compares rates")
    mi2_case(multigrid           contract "grid levels: the refusals and MOSE's V cycle on the finest grid")
    mi2_case(chimera             contract "an overset mesh: one warning with MOSE's chimera-cell count")
    mi2_case(ice-bc              exact    "ICE boundary types 301 and 403 and the population count")
    mi2_case(mass-flux-inlet     exact    "the gas inflow handed to IGLOO at 401 and 403 inlets")
    mi2_case(restart             identity "MOSE + ICE restarts equal the run in one go; the stops and their way out")
    # fpe builds its binary in a tree that links bin/hydra_MI2 too (restored after): nothing runs beside it
    mi2_case(fpe                 robustness "igloo-evap at gas-order 2 under -fpe0 (builds its own binary)" SERIAL)

    foreach(_v check_jpl evap-box/check JPL-3solver-hot-2micron/check)
        string(REGEX REPLACE "/check$" "" _n "${_v}")
        string(REPLACE "check_" "" _n "${_n}")
        set(_d "the JPL nozzle against Chang et al.")
        if(_n STREQUAL "evap-box")
            set(_d "the evaporating boxes against an independent integration of the drop equations")
        elseif(_n STREQUAL "JPL-3solver-hot-2micron")
            set(_d "alumina freezing in the hot nozzle: ICE against IGLOO")
        endif()
        hydra_test(MI2.validation.${_n} KIND reference DESC "${_d}" LABELS MI2 validation
            COMMAND bash "${_val}/${_v}.sh" "${_runs}/validation/${_n}" "${_exe}")
    endforeach()
endif()

# ---------------------------------------------------------------------------
#  Every line of the tiers file of this application must name a registered test
# ---------------------------------------------------------------------------
get_property(_hydra_names GLOBAL PROPERTY HYDRA_TEST_NAMES)
foreach(_n IN LISTS _tier_names)
    if(_n MATCHES "^${EXE_NAME}\\." AND NOT _n IN_LIST _hydra_names)
        message(FATAL_ERROR "test/tiers.txt lists '${_n}', which no hydra_test() registers")
    endif()
endforeach()
