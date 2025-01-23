﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ODMR_Lab.IO操作;

namespace ODMR_Lab.实验部分.场效应器件测量
{
    public class IVMeasureExpParams : ExpParamBase
    {
        public Param<string> DeviceName { get; set; } = new Param<string>("设备名称", "");
    }
}
