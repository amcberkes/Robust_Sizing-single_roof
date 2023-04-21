## Sample JSON object and specification

```
{
  "MyLocation": {
    "city": "string",
    "latitude": "number",
    "longitude": "number"
  },
  "SolarPanelParameters": {
    "systemLosses": "number" // between 0 - 99
  },
  "SolarPanelPosition": {
    "numberOfRoofSegments": "number",
    "roofSegments": [
      {
        "tilt": "number", // between -90 to 90, accurate to 2 decimal places
        "azimuth": "number" //between 0 to 180, accurate to 2 decimal places
      }
    ]
  },
  "SolarPanelType": {
    "moduleType": "string", // one of 'Standard', 'Premium', 'Thin film'
    "arrayType": "string" // one of 'Fixed - Open Rack', 'Fixed - Roof Mounted', '1 - Axis', '1 - Axis Backtracking', '2 - Axis'
  },
  "ElectricityLoadEstimation": {
    "numVehicles": "number",
    "bidirectionalVehicles": "number",
    "batteryCapacity": [
      {
        "vehicle": "number",
        "capacity": "number" // battery capacity for each vehicle in kWh
      }
    ],
    "stateOfCharge": "number", // "When start recharging, what is the state of charge you usually re-charge to?" - between 0 to 100
    "monthlyElectricityLoad": [
      {
        "month": "number",
        "load": "number"
      }
    ],
    "weeklyCommutingTable": {
      "monday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string", // in the format of HH:mm
          "returnAt": "string",  // in the format of HH:mm
          "distance": "number",
          "stateOfCharge": "number"
        }
      ],
      "tuesday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string",
          "returnAt": "string",
          "distance": "number",
          "stateOfCharge": "number"
        }
      ],
      "wednesday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string",
          "returnAt": "string",
          "distance": "number",
          "stateOfCharge": "number"
        }
      ],
      "thursday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string",
          "returnAt": "string",
          "distance": "number",
          "stateOfCharge": "number"
        }
      ],
      "friday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string",
          "returnAt": "string",
          "distance": "number",
          "stateOfCharge": "number"
        }
      ],
      "saturday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string",
          "returnAt": "string",
          "distance": "number",
          "stateOfCharge": "number"
        }
      ],
      "sunday": [
        {
          "vehicleNo": "number",
          "leaveAt": "string",
          "returnAt": "string",
          "distance": "number",
          "stateOfCharge": "number"
        }
      ]
    }
  },
  "costs": {
    "solarPanelPrice": "number", // it will be a number between 500 to 5000
    "batteryPrice": "number" // it will be a number between 0 to 2500
  },
  "capacities": {
    "solarPanelCapacity": "number", // it will be a number greater than 1 and default to 10
    "batteryCapacity": "number" // it will be a number greater than 1 and default to 20
  },
  "confidenceLevel": "number", // it will be a number between 0.5 to 1
  "daysInSample": "number" // it will be a number between 1 to 365
}
```