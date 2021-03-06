# RESA
Repeating Event Sequence Alarm

The RESA was designed to alert users to sequences of repeating events that sometimes occur during volcanic activity. These sequences may precede eruptions, so the RESA was created to provide automatic early detection and short-term warning of possible eruptions. However, the RESA is capable of detecting and alerting on any type of repeating activity, including that sometimes preceding landslides/avalanches and repetitive small explosions.

The idea behind the repeating event sequence alarm (RESA) is to use a standard detector (e.g., STA/LTA) to find events and then apply a correlation-matching algorithm to identify repeating event sequences. Once a sequence has been identified in progress on a minimum number of stations, an alert notification (email or text) is sent to specified recipients. This procedure can also be followed for a “level 2” alert that notifies of an increase in event rate. Notifications are also sent when a sequence event rate drops below the minimum requirements and is deemed to have ended. The alarm runs periodically on a chosen time interval.

To reduce the computational burden of the RESA, only a specified maximum number of template waveforms are saved for each station at any given time. The low computational requirements are intended to create an alarm that can be run in the background for long periods of time with little interference in the computer’s day-to-day processes. 

Visit the Wiki for more details.

For more information (and citation!), you can also read the paper here: Tepp, G. (2018). A Repeating Event Sequence Alarm for Monitoring Volcanoes. Seismological Research Letters, doi:10.1785/0220170263.

Note: The RESA Matlab code uses GISMO (available here: https://github.com/geoscience-community-codes/GISMO) with a couple edits to files available here. There are also plans to make a Python version of at least the real-time RESA.
