﻿using CodeHelper;
using HardWares.纳米位移台;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ODMR_Lab.位移台部分
{

    /// <summary>
    /// 位移台信息，包含位移台方向和是否反向的信息
    /// </summary>
    public class NanoMoverInfo : DeviceInfoBase<NanoControllerBase>
    {
        public List<NanoStageInfo> Stages { get; private set; } = new List<NanoStageInfo>();

        public NanoMoverInfo()
        {
        }

        public override DeviceInfoBase<NanoControllerBase> CreateDeviceInfo(NanoControllerBase device, DeviceConnectInfo connectinfo)
        {
            NanoMoverInfo info = new NanoMoverInfo() { Device = device, ConnectInfo = connectinfo };
            foreach (var item in info.Device.Stages)
            {
                info.Stages.Add(new NanoStageInfo(info, item));
            }
            return info;
        }



        protected override void AutoConnectedAction(FileObject file)
        {
            foreach (var item in Stages)
            {
                if (file.Descriptions.Keys.Contains("AxisPart " + item.Device.AxisName))
                {
                    item.PartType = (PartTypes)Enum.Parse(typeof(PartTypes), file.Descriptions["AxisPart " + item.Device.AxisName]);
                }
                if (file.Descriptions.Keys.Contains("AxisType " + item.Device.AxisName))
                {
                    item.MoverType = (MoverTypes)Enum.Parse(typeof(MoverTypes), file.Descriptions["AxisType " + item.Device.AxisName]);
                }
            }
        }

        protected override void CloseFileAction(FileObject obj)
        {
            foreach (var item in Stages)
            {
                obj.Descriptions.Add("AxisPart " + item.Device.AxisName, Enum.GetName(item.PartType.GetType(), item.PartType));
                obj.Descriptions.Add("AxisType " + item.Device.AxisName, Enum.GetName(item.MoverType.GetType(), item.MoverType));
            }
        }
    }

    public class NanoStageInfo
    {
        public MoverTypes MoverType { get; set; } = MoverTypes.None;

        /// <summary>
        /// 所属部分
        /// </summary>
        public PartTypes PartType { get; set; } = PartTypes.None;

        public NanoMoverInfo Parent { get; set; } = null;

        public StageBase Device { get; set; } = null;

        public NanoStageInfo(NanoMoverInfo info, StageBase stage)
        {
            Parent = info;
            Device = stage;
        }

        private object lockobj = new object();

        public bool IsWriting { get; private set; } = false;

        public bool Use(bool showmessagebox = false, bool log = true)
        {
            lock (lockobj)
            {
                if (IsWriting)
                {
                    MessageLogger.AddLogger("设备", "未能成功获取位移台设备" + Parent.Device.ProductName + "的" + Device.AxisName + "轴,轴正在使用。", MessageTypes.Warning, showmessagebox, log);
                    return false;
                }
                IsWriting = true;
                return true;
            }
        }

        public void UnUse()
        {
            lock (lockobj)
                IsWriting = false;
        }
    }


    /// <summary>
    /// 移动类型
    /// </summary>
    public enum MoverTypes
    {
        /// <summary>
        /// 未分配位移台
        /// </summary>
        None = -1,
        /// <summary>
        /// X位移台
        /// </summary>
        X = 0,
        /// <summary>
        /// Y位移台
        /// </summary>
        Y = 1,
        /// <summary>
        /// Z位移台
        /// </summary>
        Z = 2,
        /// <summary>
        /// X角度偏摆台
        /// </summary>
        AngleX = 3,
        /// <summary>
        /// Y角度偏摆台
        /// </summary>
        AngleY = 4,
        /// <summary>
        /// Z角度偏摆台
        /// </summary>
        AngleZ = 5
    }

    /// <summary>
    /// 所属部分类型
    /// </summary>
    public enum PartTypes
    {
        /// <summary>
        /// 未分配
        /// </summary>
        None = -1,
        /// <summary>
        /// 探针部分
        /// </summary>
        Probe = 0,
        /// <summary>
        /// 样品部分
        /// </summary>
        Sample = 1,
        /// <summary>
        /// 微波部分
        /// </summary>
        Microwave = 2,
        /// <summary>
        /// 磁铁部分
        /// </summary>
        Magnnet = 3,
    }
}
