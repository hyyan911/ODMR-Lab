﻿using CodeHelper;
using HardWares.仪器列表.电动翻转座;

namespace ODMR_Lab.设备部分.相机
{
    public class LockinInfo : DeviceInfoBase<FlipMotorBase>
    {
        public override void CreateDeviceInfoBehaviour()
        {
        }

        public override string GetDeviceDescription()
        {
            return Device.ProductName;
        }

        protected override void AutoConnectedAction(FileObject file)
        {
        }

        protected override void CloseFileAction(FileObject obj)
        {
        }
    }
}
